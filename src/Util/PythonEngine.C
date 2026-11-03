/*******************************************************************************

  Copyright (C) 2026 Roberto Peverati
  Developed with the assistance of Claude (Anthropic).

  This file is part of IQmol, a free molecular visualization program. See
  <http://iqmol.org> for more details.

  IQmol is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  IQmol is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along
  with IQmol.  If not, see <http://www.gnu.org/licenses/>.

********************************************************************************/

#include "PythonEngine.h"
#include "Optimizer.h"
#include "Vibrations.h"
#include <QCoreApplication>
#include <QDir>
#include <QElapsedTimer>
#include <QJsonDocument>
#include <QProcess>
#include <QProcessEnvironment>
#include <cmath>


namespace IQmol {
namespace Calc {

namespace {

double const s_bohrInAngstrom(0.529177210903);
double const s_eVPerHartree(27.211386245988);
double const s_debyePerEAngstrom(4.80320471);

// The forces of the helper (eV/Angstrom) are the negative gradient in Hartree/Bohr
double const s_gradientFactor(s_bohrInAngstrom / s_eVPerHartree);
double const s_hessianFactor(s_bohrInAngstrom * s_bohrInAngstrom / s_eVPerHartree);

} // end anonymous namespace


/// Why a request to the helper did not give a result
struct PythonFailure {
   QString code;        // error code of the protocol, "cancelled", or "helper" for a lost helper
   QString message;
   QString detail;

   bool cancelled() const { return code == "cancelled"; }
   bool failed() const { return !code.isEmpty(); }
   void set(QString const& c, QString const& m, QString const& d = QString()) 
   { 
      code = c; message = m; detail = d; 
   }
};


/// The helper process and the JSON lines protocol, see python/iqmolx-engines/PROTOCOL.md.
/// Requests are synchronous: request() returns when the helper has answered.
class PythonHelper {
   public:
      PythonHelper() : m_process(0), m_nextId(1) { }
      ~PythonHelper() { stop(); }

      bool running() const { return m_process && m_process->state() == QProcess::Running; }

      bool start(QString const& python, QString const& pythonPath, bool interactive, 
         PythonFailure& failure)
      {
         stop();
         m_stderr.clear();
         m_process = new QProcess();
         m_process->setProcessChannelMode(QProcess::SeparateChannels);

         QProcessEnvironment env(QProcessEnvironment::systemEnvironment());
         if (!pythonPath.isEmpty()) {
            QString old(env.value("PYTHONPATH"));
            env.insert("PYTHONPATH", old.isEmpty() ? pythonPath 
               : pythonPath + QDir::listSeparator() + old);
         }
         m_process->setProcessEnvironment(env);

         QString program(python.isEmpty() ? QString("python3") : python);
         m_process->start(program, QStringList() << "-u" << "-m" << "iqmolx_engines");
         if (!m_process->waitForStarted(10000)) {
            failure.set("helper", "Unable to start the Python interpreter '" + program + "'",
               m_process->errorString());
            stop();
            return false;
         }

         // The helper announces itself, then both sides agree on the protocol version
         QElapsedTimer clock;
         clock.start();
         while (true) {
            QJsonObject message;
            Read read(readMessage(message, failure, interactive));
            if (read == Lost) { stop(); return false; }
            if (read == Nothing) {
               if (clock.elapsed() > 60000) {
                  failure.set("helper", "The Python helper did not start within a minute", 
                     stderrTail());
                  stop();
                  return false;
               }
               continue;
            }
            if (message.value("event").toString() == "ready") break;
            if (message.value("event").toString() == "fatal") {
               failure.set("helper", "The Python helper cannot start: " 
                  + message.value("message").toString(), stderrTail());
               stop();
               return false;
            }
         }

         QJsonObject fields, result;
         fields["protocol"] = 1;
         fields["client"] = "IQMolX";
         if (!request("hello", fields, result, failure, 0, 0, interactive)) { stop(); return false; }
         return true;
      }

      void stop()
      {
         if (!m_process) return;
         if (m_process->state() == QProcess::Running) {
            m_process->closeWriteChannel();            // end of input ends the helper
            if (!m_process->waitForFinished(2000)) {
               m_process->kill();
               m_process->waitForFinished(2000);
            }
         }
         delete m_process;
         m_process = 0;
      }

      typedef std::function<void(double fraction, QString const& message)> Progress;

      /// Sends a request and waits for its final message.  Returns false and sets 'failure'
      /// when the helper answers with an error, is lost, or *cancel becomes true; in the last
      /// case the helper is told to cancel and the answer is no longer awaited.
      bool request(QString const& type, QJsonObject fields, QJsonObject& result, 
         PythonFailure& failure, Progress const& progress, bool const* cancel, bool interactive)
      {
         failure = PythonFailure();
         if (!running()) {
            failure.set("helper", "The Python helper is not running", stderrTail());
            return false;
         }

         int id(m_nextId++);
         fields["id"] = id;
         fields["type"] = type;
         send(fields);

         while (true) {
            if (cancel && *cancel) {
               QJsonObject message;
               message["id"] = m_nextId++;
               message["type"] = "cancel";
               message["target"] = id;
               send(message);
               failure.set("cancelled", "Calculation cancelled");
               return false;
            }

            QJsonObject message;
            Read read(readMessage(message, failure, interactive));
            if (read == Lost) return false;
            if (read == Nothing) continue;
            if (message.value("id").toInt(-1) != id) continue;   // late news of an abandoned request

            if (message.value("event").toString() == "progress") {
               if (progress) progress(message.value("fraction").toDouble(-1.0), 
                  message.value("message").toString());
               continue;
            }
            if (message.value("ok").toBool()) {
               result = message.value("result").toObject();
               return true;
            }
            QJsonObject error(message.value("error").toObject());
            failure.set(error.value("code").toString("engine_error"), 
               error.value("message").toString("The Python helper reported an error"),
               error.value("detail").toString());
            return false;
         }
      }

      QString stderrTail() const { return m_stderr.right(3000).trimmed(); }

   private:
      void send(QJsonObject const& message)
      {
         m_process->write(QJsonDocument(message).toJson(QJsonDocument::Compact) + "\n");
      }

      void collectStderr()
      {
         m_stderr += QString::fromLocal8Bit(m_process->readAllStandardError());
         if (m_stderr.size() > 20000) m_stderr = m_stderr.right(10000);
      }

      enum Read { Got, Nothing, Lost };

      /// Looks for the next complete line of the helper as an object.  Waits a short while
      /// if there is none (Nothing), and gives Lost if the helper is gone.
      Read readMessage(QJsonObject& message, PythonFailure& failure, bool interactive)
      {
         collectStderr();
         while (m_process->canReadLine()) {
            QJsonParseError parseError;
            QJsonDocument doc(QJsonDocument::fromJson(m_process->readLine(), &parseError));
            if (parseError.error == QJsonParseError::NoError && doc.isObject()) {
               message = doc.object();
               return Got;
            }
         }
         if (m_process->state() == QProcess::NotRunning) {
            m_process->waitForReadyRead(0);
            collectStderr();
            if (m_process->canReadLine()) return readMessage(message, failure, interactive);
            QString how(m_process->exitStatus() == QProcess::CrashExit ? "crashed" 
               : "stopped (exit code " + QString::number(m_process->exitCode()) + ")");
            failure.set("helper", "The Python helper " + how, stderrTail());
            return Lost;
         }
         m_process->waitForReadyRead(20);
         if (interactive) QCoreApplication::processEvents();
         return Nothing;
      }

      QProcess* m_process;
      int       m_nextId;
      QString   m_stderr;
};


PythonEngine::PythonEngine(QString const& python, QString const& engine) 
 : m_python(python), m_engine(engine), m_helper(new PythonHelper()), 
   m_analyticHessian(true), m_hessianStep(0.01), m_interactive(false), m_cancelled(false)
{
}


PythonEngine::~PythonEngine()
{
   delete m_helper;
}


QString PythonEngine::name() const
{
   QString title(m_descriptor.value("title").toString());
   if (!title.isEmpty()) return title;
   return m_engine.isEmpty() ? QString("Python engine") : m_engine;
}


void PythonEngine::setPython(QString const& path)
{
   if (path == m_python) return;
   m_python = path;
   shutdown();
}


void PythonEngine::setPythonPath(QString const& path)
{
   if (path == m_pythonPath) return;
   m_pythonPath = path;
   shutdown();
}


void PythonEngine::setEngine(QString const& name)
{
   if (name == m_engine) return;
   m_engine = name;
   m_descriptor = QJsonObject();
}


void PythonEngine::setLoadOptions(QJsonObject const& options) { m_loadOptions = options; }
void PythonEngine::setComputeOptions(QJsonObject const& options) { m_computeOptions = options; }


void PythonEngine::shutdown()
{
   m_helper->stop();
   m_loadedKey.clear();
   m_engines = QJsonArray();
   m_descriptor = QJsonObject();
}


void PythonEngine::report(QString const& text)
{
   m_log += text;
   if (m_progress) m_progress(text);
}


QJsonArray PythonEngine::listEngines(QString* error)
{
   PythonFailure failure;
   bool interactive(m_interactive);
   if (!m_helper->running()) {
      m_loadedKey.clear();
      if (!m_helper->start(m_python, m_pythonPath, interactive, failure)) {
         if (error) *error = failure.message + (failure.detail.isEmpty() ? "" : "\n" + failure.detail);
         return QJsonArray();
      }
   }

   QJsonObject result;
   if (!m_helper->request("list_engines", QJsonObject(), result, failure, 0, 0, interactive)) {
      if (error) *error = failure.message;
      return QJsonArray();
   }
   m_engines = result.value("engines").toArray();
   for (QJsonValue const& value : m_engines) {
      if (value.toObject().value("name").toString() == m_engine) m_descriptor = value.toObject();
   }
   return m_engines;
}


bool PythonEngine::offers(QString const& property) const
{
   for (QJsonValue const& value : m_descriptor.value("properties").toArray()) {
      if (value.toString() == property) return true;
   }
   return false;
}


/// Makes sure the helper runs and holds the model that was asked for
bool PythonEngine::prepare(PythonFailure& failure)
{
   if (m_engine.isEmpty()) {
      failure.set("bad_request", "No Python engine has been selected");
      return false;
   }

   if (!m_helper->running()) {
      m_loadedKey.clear();
      m_engines = QJsonArray();
      if (!m_helper->start(m_python, m_pythonPath, m_interactive, failure)) return false;
   }

   if (m_engines.isEmpty()) {
      QString error;
      listEngines(&error);
      if (m_engines.isEmpty()) {
         failure.set("helper", "Unable to list the engines of the Python helper\n" + error);
         return false;
      }
   }
   if (m_descriptor.isEmpty()) {
      failure.set("unknown_engine", "The Python helper has no engine '" + m_engine + "'");
      return false;
   }

   QString key(m_engine + QChar('\n') + 
      QString::fromUtf8(QJsonDocument(m_loadOptions).toJson(QJsonDocument::Compact)));
   if (key == m_loadedKey) return true;

   report("Loading " + name() + " ...\n");
   QJsonObject fields, result;
   fields["engine"] = m_engine;
   fields["options"] = m_loadOptions;
   auto progress = [this](double fraction, QString const& message) {
      if (fraction >= 0.0) report(QString::asprintf("  %3.0f%%  ", 100.0*fraction) + message + "\n");
      else report("  " + message + "\n");
   };
   if (!m_helper->request("load", fields, result, failure, progress, &m_cancelled, m_interactive)) {
      return false;
   }

   QJsonObject info(result.value("info").toObject());
   QStringList parts;
   for (QString const& k : info.keys()) parts << k + ": " + info.value(k).toVariant().toString();
   if (!parts.isEmpty()) report("  " + parts.join(", ") + "\n");
   m_loadedKey = key;
   return true;
}


/// One compute request.  The coordinates are in Bohr.
bool PythonEngine::compute(Input const& input, std::vector<double> const& bohr, 
   QStringList const& want, QJsonObject& result, PythonFailure& failure)
{
   QJsonArray numbers, coordinates, properties;
   for (int z : input.atomicNumbers) numbers.append(z);
   for (size_t i = 0; i + 2 < bohr.size(); i += 3) {
      QJsonArray atom;
      for (int k = 0; k < 3; ++k) atom.append(bohr[i+k] * s_bohrInAngstrom);
      coordinates.append(atom);
   }
   for (QString const& w : want) properties.append(w);

   QJsonObject fields;
   fields["numbers"] = numbers;
   fields["coordinates"] = coordinates;
   fields["charge"] = input.charge;
   fields["multiplicity"] = input.multiplicity;
   fields["want"] = properties;
   fields["options"] = m_computeOptions;
   return m_helper->request("compute", fields, result, failure, 0, &m_cancelled, m_interactive);
}


/// Energy (Hartree) and gradient (Hartree/Bohr) at a geometry in Bohr
bool PythonEngine::gradient(Input const& input, std::vector<double> const& bohr, double& energy,
   std::vector<double>& gradient, PythonFailure& failure, QJsonObject* details)
{
   QStringList want;
   want << "energy" << "forces";
   if (details && offers("dipole")) want << "dipole";

   QJsonObject result;
   if (!compute(input, bohr, want, result, failure)) return false;

   energy = result.value("energy").toDouble() / s_eVPerHartree;
   gradient.assign(bohr.size(), 0.0);
   QJsonArray forces(result.value("forces").toArray());
   for (int i = 0; i < forces.size(); ++i) {
      QJsonArray f(forces[i].toArray());
      for (int k = 0; k < 3 && k < f.size(); ++k) {
         gradient[3*i+k] = -f[k].toDouble() * s_gradientFactor;
      }
   }
   if (details) *details = result;
   return true;
}


/// Hessian in Hartree/Bohr^2, row major: the one of the engine if it has one, otherwise by 
/// central differences of the forces
bool PythonEngine::hessian(Input const& input, std::vector<double> const& bohr, 
   std::vector<double>& hessian, PythonFailure& failure)
{
   size_t const n(bohr.size());
   hessian.assign(n*n, 0.0);

   if (m_analyticHessian && offers("hessian")) {
      report("Hessian from " + name() + "\n");
      QJsonObject result;
      if (!compute(input, bohr, QStringList("hessian"), result, failure)) return false;
      QJsonArray rows(result.value("hessian").toArray());
      for (int i = 0; i < rows.size() && i < static_cast<int>(n); ++i) {
         QJsonArray row(rows[i].toArray());
         for (int j = 0; j < row.size() && j < static_cast<int>(n); ++j) {
            hessian[i*n+j] = row[j].toDouble() * s_hessianFactor;
         }
      }
      for (size_t i = 0; i < n; ++i) {
         for (size_t j = i+1; j < n; ++j) {
            hessian[i*n+j] = hessian[j*n+i] = 0.5 * (hessian[i*n+j] + hessian[j*n+i]);
         }
      }
      return true;
   }

   report(QString::asprintf("Numerical Hessian from forces (%d displaced calculations) ...\n", 
      static_cast<int>(2*n)));
   double const h(m_hessianStep);
   for (size_t i = 0; i < n; ++i) {
      std::vector<double> x(bohr), plus, minus;
      double e;
      x[i] = bohr[i] + h;
      if (!gradient(input, x, e, plus, failure)) return false;
      x[i] = bohr[i] - h;
      if (!gradient(input, x, e, minus, failure)) return false;
      for (size_t j = 0; j < n; ++j) hessian[i*n+j] = (plus[j] - minus[j]) / (2.0*h);
   }
   for (size_t i = 0; i < n; ++i) {
      for (size_t j = i+1; j < n; ++j) {
         hessian[i*n+j] = hessian[j*n+i] = 0.5 * (hessian[i*n+j] + hessian[j*n+i]);
      }
   }
   return true;
}


bool PythonEngine::validInput(Input const& input, Result& result)
{
   int natoms(input.atomicNumbers.size());
   if (natoms == 0 || input.coordinates.size() != 3*natoms) {
      result.message = "Invalid input: inconsistent number of atoms and coordinates";
      return false;
   }
   if (input.multiplicity < 1) {
      result.message = "Invalid spin multiplicity";
      return false;
   }
   return true;
}


Result PythonEngine::failed(PythonFailure const& failure, QString const& log)
{
   Result result;
   result.log = log;
   if (failure.cancelled()) {
      result.cancelled = true;
      result.message = "Calculation cancelled";
      return result;
   }
   result.message = failure.message;
   if (!failure.detail.isEmpty()) result.message += "\n" + failure.detail;
   return result;
}


static std::vector<double> toBohr(QVector<double> const& angstrom)
{
   std::vector<double> bohr(angstrom.size());
   for (int i = 0; i < angstrom.size(); ++i) bohr[i] = angstrom[i] / s_bohrInAngstrom;
   return bohr;
}


static double norm(std::vector<double> const& v)
{
   double sum(0.0);
   for (double x : v) sum += x*x;
   return std::sqrt(sum);
}


static double dipoleDebye(QJsonObject const& details)
{
   QJsonArray d(details.value("dipole").toArray());
   if (d.size() != 3) return -1.0;
   return s_debyePerEAngstrom * std::sqrt(d[0].toDouble()*d[0].toDouble() 
        + d[1].toDouble()*d[1].toDouble() + d[2].toDouble()*d[2].toDouble());
}


Result PythonEngine::energy(Input const& input)
{
   Result result;
   m_cancelled = false;
   m_log.clear();
   if (!validInput(input, result)) return result;

   QElapsedTimer clock;
   clock.start();

   PythonFailure failure;
   if (!prepare(failure)) return failed(failure, m_log);

   std::vector<double> x(toBohr(input.coordinates)), g;
   double e(0.0);
   QJsonObject details;
   if (!gradient(input, x, e, g, failure, &details)) return failed(failure, m_log);

   result.energy = e;
   result.initialEnergy = e;
   result.gradientNorm = norm(g);
   result.dipoleDebye = dipoleDebye(details);
   result.coordinates = input.coordinates;
   report(QString::asprintf("%s\n  energy %.8f Eh   |g| = %.6f Eh/Bohr\n", 
      qPrintable(name()), e, result.gradientNorm));

   result.converged = true;
   result.seconds = clock.elapsed() / 1000.0;
   result.log = m_log;
   result.success = true;
   return result;
}


Result PythonEngine::optimize(Input const& input)
{
   Result result;
   m_cancelled = false;
   m_log.clear();
   if (!validInput(input, result)) return result;

   QElapsedTimer clock;
   clock.start();

   PythonFailure failure;
   if (!prepare(failure)) return failed(failure, m_log);

   report("Geometry optimization with " + name() + " and IQmol's L-BFGS optimizer\n"
          "convergence: energy change < 5e-6 Eh, gradient norm < 1e-3 Eh/Bohr\n\n");

   bool first(true);
   auto evaluate = [&](std::vector<double> const& xx, double& e, std::vector<double>& g) {
      bool ok(gradient(input, xx, e, g, failure));
      if (ok && first) { result.initialEnergy = e; first = false; }
      return ok;
   };

   auto progress = [&](int cycle, double e, double gnorm, std::vector<double> const&) {
      report(QString::asprintf("  cycle %3d   E = %.8f Eh   |g| = %.6f Eh/Bohr\n", cycle, e, gnorm));
      if (m_interactive) QCoreApplication::processEvents();
      return !m_cancelled;
   };

   std::vector<double> x(toBohr(input.coordinates));
   OptimizerResult opt(Optimizer().minimize(evaluate, x, progress));
   result.cycles = opt.cycles;
   result.converged = opt.converged;

   if (opt.cancelled || failure.cancelled()) {
      result.cancelled = true;
      result.message = "Calculation cancelled";
      result.log = m_log;
      return result;
   }
   if (opt.x.empty() || (!opt.converged && failure.failed())) {
      Result error(failed(failure, m_log));
      error.cycles = result.cycles;
      return error;
   }

   x = opt.x;
   report(opt.converged ? "\n*** GEOMETRY OPTIMIZATION CONVERGED AFTER " 
                        + QString::number(opt.cycles) + " CYCLES ***\n\n"
                      : "\n*** " + QString::fromStdString(opt.message) + " ***\n\n");
   if (!opt.converged) result.message = "Geometry optimization did not converge\n" 
                                       + QString::fromStdString(opt.message);

   // A last calculation at the result for the details
   std::vector<double> g;
   double e(0.0);
   QJsonObject details;
   if (!gradient(input, x, e, g, failure, &details)) {
      Result error(failed(failure, m_log));
      error.cycles = result.cycles;
      return error;
   }
   result.energy = e;
   result.gradientNorm = norm(g);
   result.dipoleDebye = dipoleDebye(details);
   report(QString::asprintf("  final energy %.8f Eh   |g| = %.6f Eh/Bohr\n", e, result.gradientNorm));

   result.coordinates.resize(x.size());
   for (size_t i = 0; i < x.size(); ++i) result.coordinates[i] = x[i] * s_bohrInAngstrom;

   result.seconds = clock.elapsed() / 1000.0;
   result.log = m_log;
   result.success = true;
   return result;
}


Result PythonEngine::frequencies(Input const& input)
{
   Result result;
   m_cancelled = false;
   m_log.clear();
   if (!validInput(input, result)) return result;

   QElapsedTimer clock;
   clock.start();

   PythonFailure failure;
   if (!prepare(failure)) return failed(failure, m_log);

   // The energy and gradient of the geometry: the frequencies only mean something at a 
   // stationary point
   std::vector<double> x(toBohr(input.coordinates)), g;
   double e(0.0);
   QJsonObject details;
   if (!gradient(input, x, e, g, failure, &details)) return failed(failure, m_log);
   result.energy = e;
   result.initialEnergy = e;
   result.gradientNorm = norm(g);
   result.dipoleDebye = dipoleDebye(details);
   report(QString::asprintf("%s\n  energy %.8f Eh   |g| = %.6f Eh/Bohr\n\n", 
      qPrintable(name()), e, result.gradientNorm));

   std::vector<double> h;
   if (!hessian(input, x, h, failure)) return failed(failure, m_log);

   std::vector<int> numbers(input.atomicNumbers.toStdVector());
   VibrationalAnalysis analysis(analyseVibrations(numbers, x, h));
   if (!analysis.ok) {
      result.message = "Vibrational analysis failed: " + QString::fromStdString(analysis.message);
      result.log = m_log;
      return result;
   }

   VibrationalData& data(result.vibrations);
   for (size_t i = 0; i < analysis.frequencies.size(); ++i) {
       data.frequencies.append(analysis.frequencies[i]);
       data.reducedMasses.append(analysis.reducedMasses[i]);
       data.modes.append(QVector<double>::fromStdVector(analysis.modes[i]));
   }
   data.zeroPointEnergy = analysis.zeroPointEnergy;
   result.hasVibrations = true;

   report("\n  mode     frequency (cm-1)\n");
   for (int i = 0; i < data.frequencies.size(); ++i) {
       report(QString::asprintf("  %4d     %12.2f\n", i + 1, data.frequencies[i]));
   }
   report(QString::asprintf("\n  zero point energy  %.8f Eh\n", data.zeroPointEnergy));

   result.coordinates = input.coordinates;
   result.converged = true;
   result.seconds = clock.elapsed() / 1000.0;
   result.log = m_log;
   result.success = true;
   return result;
}

} } // end namespace IQmol::Calc
