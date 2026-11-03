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

#include "XtbProcess.h"
#include "XtbOutput.h"
#include <QCoreApplication>
#include <QDir>
#include <QElapsedTimer>
#include <QEventLoop>
#include <QTimer>
#include <QFile>
#include <QFileInfo>
#include <QProcess>
#include <QProcessEnvironment>
#include <QRegularExpression>
#include <QStandardPaths>
#include <QTemporaryDir>
#include <QTextStream>


namespace IQmol {
namespace Calc {

namespace {

// Element symbols for atomic numbers 0-118
char const* s_symbols[] = {
   "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg",
   "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn",
   "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
   "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
   "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
   "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir",
   "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
   "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
   "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv",
   "Ts", "Og"
};

int const s_numSymbols(sizeof(s_symbols)/sizeof(s_symbols[0]));


QString xyzText(Input const& input)
{
   QString text;
   QTextStream out(&text);
   out << input.atomicNumbers.size() << "\n\n";
   for (int i = 0; i < input.atomicNumbers.size(); ++i) {
       int Z(input.atomicNumbers[i]);
       out << s_symbols[(Z > 0 && Z < s_numSymbols) ? Z : 0] << " "
           << QString::number(input.coordinates[3*i+0], 'f', 10) << " "
           << QString::number(input.coordinates[3*i+1], 'f', 10) << " "
           << QString::number(input.coordinates[3*i+2], 'f', 10) << "\n";
   }
   return text;
}


// Reads the Cartesian coordinates from an xyz file.  Recent versions of xtb write
// extended xyz with a long comment line, but the atom lines are standard.
bool readXyz(QString const& path, int natoms, QVector<double>& coordinates)
{
   QFile file(path);
   if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return false;

   QTextStream in(&file);
   bool ok(false);
   int n(in.readLine().trimmed().toInt(&ok));
   if (!ok || n != natoms) return false;
   in.readLine();  // comment

   coordinates.clear();
   for (int i = 0; i < natoms; ++i) {
       QStringList tokens(in.readLine().split(QRegularExpression("\\s+"),
          Qt::SkipEmptyParts));
       if (tokens.size() < 4) return false;
       for (int j = 1; j < 4; ++j) {
           double x(tokens[j].toDouble(&ok));
           if (!ok) return false;
           coordinates.append(x);
       }
   }
   return true;
}


QString tail(QString const& text, int lines = 15)
{
   QStringList list(text.split('\n'));
   while (!list.isEmpty() && list.last().trimmed().isEmpty()) list.removeLast();
   return list.mid(qMax(0, list.size()-lines)).join('\n');
}


// xtb uses the first of a repeated option, so an option given by the user has to replace the
// one IQmol sets rather than follow it.  Only the options that the xtb window itself sets
// (method, charge, spin, solvent, electronic temperature) are replaced; the kind of job
// is left alone, as the results are read according to it.
void applyUserOptions(QStringList& args, QStringList const& user)
{
   struct Option { QStringList names; bool hasValue; };
   static QList<Option> const options = {
      { { "--gfn", "--gfn0", "--gfn1", "--gfn2" }, true },  // "--gfn 2" and "--gfn2" 
      { { "--gfnff" },                             false },
      { { "--chrg", "-c" },                        true },
      { { "--uhf", "-u" },                         true },
      { { "--alpb", "--gbsa" },                    true },
      { { "--etemp" },                             true }
   };

   // the method: "--gfn 2", "--gfn2" and "--gfnff" all select one
   auto isMethod = [](QString const& name) { return name.startsWith("--gfn"); };

   for (QString const& token : user) {
       if (!token.startsWith('-')) continue;
       for (Option const& option : options) {
           if (!option.names.contains(token)) continue;
           // remove ours
           for (Option const& ours : options) {
               bool same(isMethod(token) ? isMethod(ours.names.first()) 
                                         : ours.names.contains(token) || 
                                           (token == "--gbsa" && ours.names.contains("--alpb")));
               if (!same) continue;
               for (QString const& name : ours.names) {
                   int i(args.indexOf(name));
                   if (i >= 0) args.removeAt(i);
                   if (i >= 0 && ours.hasValue && i < args.size()) args.removeAt(i);
               }
           }
       }
   }
   args << user;
}

} // end anonymous namespace


XtbProcess::XtbProcess(QString const& executable, Method method)
 : m_executable(executable), m_method(method), m_timeout(600), m_interactive(false),
   m_keepTemporary(false),
   m_cancelled(false), m_process(0)
{
}


void XtbProcess::cancel()
{
   m_cancelled = true;
   if (m_process) m_process->kill();
}


QString XtbProcess::name() const
{
   switch (m_method) {
      case GFN2:  return "GFN2-xTB";
      case GFN1:  return "GFN1-xTB";
      case GFNFF: return "GFN-FF";
   }
   return "xTB";
}


QStringList XtbProcess::methodArguments() const
{
   switch (m_method) {
      case GFN2:  return QStringList() << "--gfn" << "2";
      case GFN1:  return QStringList() << "--gfn" << "1";
      case GFNFF: return QStringList() << "--gfnff";
   }
   return QStringList();
}


QString XtbProcess::findExecutable(QString const& hint)
{
   QStringList candidates;
   if (!hint.isEmpty()) candidates << hint;
   QString env(qEnvironmentVariable("XTB_PATH"));
   if (!env.isEmpty()) candidates << env;

   for (QString const& path : candidates) {
       QFileInfo info(path);
       if (info.isDir()) info.setFile(QDir(path).filePath("xtb"));
       if (info.isFile() && info.isExecutable()) return info.absoluteFilePath();
   }

   // The xtb program that is delivered inside the application
#ifdef Q_OS_WIN
   QFileInfo bundled(QDir(QCoreApplication::applicationDirPath()).filePath("xtb.exe"));
#else
   QFileInfo bundled(QDir(QCoreApplication::applicationDirPath()).filePath("xtb"));
#endif
   if (bundled.isFile() && bundled.isExecutable()) return bundled.absoluteFilePath();

   // A GUI application on macOS may not inherit the shell PATH, so also
   // look in the usual locations.
   QString path(QStandardPaths::findExecutable("xtb"));
   if (path.isEmpty()) {
      path = QStandardPaths::findExecutable("xtb",
         QStringList() << "/opt/homebrew/bin" << "/usr/local/bin");
   }
   return path;
}


Result XtbProcess::energy(Input const& input)
{
   return run(input, SinglePoint);
}


Result XtbProcess::optimize(Input const& input)
{
   return run(input, Optimization);
}


Result XtbProcess::frequencies(Input const& input)
{
   return run(input, Hessian);
}


QStringList XtbProcess::arguments(Input const& input, Task task) const
{
   QStringList args;
   args << "input.xyz";
   args << methodArguments();
   args << "--chrg" << QString::number(input.charge);
   args << "--uhf"  << QString::number(input.multiplicity - 1);
   switch (task) {
      case SinglePoint:  args << "--sp";   break;
      case Optimization: args << "--opt";  break;
      case Hessian:      args << "--hess"; break;
   }
   if (!input.solvent.isEmpty()) args << "--alpb" << input.solvent;
   if (input.orbitals && m_method != GFNFF) args << "--molden";
   if (input.electronicTemperature >= 0.0 && m_method != GFNFF) {
      args << "--etemp" << QString::number(input.electronicTemperature, 'g', 10);
   }
   applyUserOptions(args, m_extraArguments);
   return args;
}


Result XtbProcess::run(Input const& input, Task task)
{
   Result result;

   int natoms(input.atomicNumbers.size());
   if (natoms == 0 || input.coordinates.size() != 3*natoms) {
      result.message = "Invalid input: inconsistent number of atoms and coordinates";
      return result;
   }
   if (input.multiplicity < 1) {
      result.message = "Invalid spin multiplicity";
      return result;
   }

   QString exe(findExecutable(m_executable));
   if (exe.isEmpty()) {
      result.message = "Unable to find the xtb executable.  Set its location in the "
                       "preferences or in the XTB_PATH environment variable.";
      return result;
   }

   // Either the directory chosen by the user, which keeps all the files, or a temporary one
   QTemporaryDir tempDir;
   QString workPath;
   QString const outName(task == Hessian ? "xtb_hess.out" : (task == Optimization ? "xtb_opt.out" : "xtb_sp.out"));
   if (m_workingDirectory.isEmpty()) {
      if (m_keepTemporary && !m_temporaryPath.isEmpty() && QDir(m_temporaryPath).exists()) {
         workPath = m_temporaryPath;     // the calculations of one job share a directory
      }else {
         if (!tempDir.isValid()) {
            result.message = "Failed to create a temporary directory";
            return result;
         }
         workPath = tempDir.path();
         if (m_keepTemporary) {
            tempDir.setAutoRemove(false);
            m_temporaryPath = workPath;
         }
      }
      if (m_keepTemporary) {
         result.directory = workPath;
         result.temporaryDirectory = true;
         // results of an earlier calculation in this directory
         QFile::remove(QDir(workPath).filePath(outName));
         if (task != Hessian) QFile::remove(QDir(workPath).filePath("xtbopt.xyz"));
         QFile::remove(QDir(workPath).filePath("molden.input"));
      }
   }else {
      workPath = QDir(m_workingDirectory).absolutePath();
      if (!QDir().mkpath(workPath)) {
         result.message = "Unable to create the working directory " + workPath;
         return result;
      }
      // Results of an earlier run must not be mistaken for those of this one
      QStringList stale;
      stale << "molden.input" << outName;
      if (task != Hessian) stale << "xtbopt.xyz" << "xtbopt.log";
      if (task == Hessian) stale << "g98.out";
      for (QString const& name : stale) {
          QFile::remove(QDir(workPath).filePath(name));
      }
      result.directory = workPath;
   }
   QDir dir(workPath);

   QString inputFile(dir.filePath("input.xyz"));
   {
      QFile file(inputFile);
      if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
         result.message = "Failed to write xtb input file";
         return result;
      }
      file.write(xyzText(input).toLocal8Bit());
   }

   QStringList args(arguments(input, task));

   QProcess process;
   process.setWorkingDirectory(workPath);
   process.setProcessChannelMode(QProcess::MergedChannels);
   QProcessEnvironment env(QProcessEnvironment::systemEnvironment());
   env.insert("OMP_STACKSIZE", "1G");
   process.setProcessEnvironment(env);

   QElapsedTimer clock;
   clock.start();
   QString output;
   bool timedOut(false);
   m_cancelled = false;

   if (m_interactive) {
      QEventLoop loop;
      QTimer timer;
      timer.setSingleShot(true);

      QObject::connect(&process, &QProcess::readyRead, [&]() {
         QString chunk(QString::fromLocal8Bit(process.readAll()));
         output += chunk;
         if (m_progress) m_progress(chunk);
      });
      QObject::connect(&process, 
         static_cast<void (QProcess::*)(int, QProcess::ExitStatus)>(&QProcess::finished),
         &loop, &QEventLoop::quit);
      QObject::connect(&process, &QProcess::errorOccurred, [&](QProcess::ProcessError error) {
         if (error == QProcess::FailedToStart) loop.quit();
      });
      QObject::connect(&timer, &QTimer::timeout, [&]() { timedOut = true; process.kill(); });

      m_process = &process;
      process.start(exe, args);
      timer.start(m_timeout * 1000);
      if (process.state() != QProcess::NotRunning || process.error() != QProcess::FailedToStart) {
         loop.exec();
      }
      m_process = 0;

      QString rest(QString::fromLocal8Bit(process.readAll()));
      output += rest;
      if (m_progress && !rest.isEmpty()) m_progress(rest);

      if (process.error() == QProcess::FailedToStart) {
         result.message = "Failed to start " + exe;
         return result;
      }
   }else {
      process.start(exe, args);

      if (!process.waitForStarted()) {
         result.message = "Failed to start " + exe;
         return result;
      }

      if (!process.waitForFinished(m_timeout * 1000)) {
         timedOut = true;
         process.kill();
         process.waitForFinished();
      }
      output = QString::fromLocal8Bit(process.readAll());
   }

   result.seconds = clock.elapsed() / 1000.0;

   if (m_cancelled) {
      result.cancelled = true;
      result.message = "Calculation cancelled";
      return result;
   }

   if (timedOut) {
      result.message = "xtb timed out after " + QString::number(m_timeout) + " seconds";
      return result;
   }

   // The complete output is kept next to the other files
   if (!m_workingDirectory.isEmpty() || m_keepTemporary) {
      QFile outFile(dir.filePath(outName));
      if (outFile.open(QIODevice::WriteOnly)) outFile.write(output.toLocal8Bit());
   }

   QFile inputXyz(inputFile);
   QString inputText;
   if (inputXyz.open(QIODevice::ReadOnly)) inputText = QString::fromLocal8Bit(inputXyz.readAll());
   result.log = exe + " " + args.join(" ") + "\n--- input.xyz ---\n" + inputText 
              + "--- xtb output ---\n" + output;

   if (process.exitStatus() != QProcess::NormalExit) {
      result.message = "xtb crashed\n" + tail(output);
      return result;
   }

   if (!parseEnergy(output, result.energy)) {
      result.message = "xtb failed (exit code " + QString::number(process.exitCode())
                     + ")\n" + tail(output);
      return result;
   }

   if (task == Optimization) {
      // xtb writes the last geometry even when the optimization does not converge
      if (!readXyz(dir.filePath("xtbopt.xyz"), natoms, result.coordinates)) {
         result.message = "Failed to read the optimized geometry from xtb\n" + tail(output);
         return result;
      }
      result.converged = output.contains("GEOMETRY OPTIMIZATION CONVERGED");
      if (!result.converged) {
         result.message = "Geometry optimization did not converge\n" + tail(output);
      }
   }else {
      result.coordinates = input.coordinates;
      result.converged = true;
   }

   parseDetails(output, result);

   if (task == Hessian) {
      // The vibrations are in the Gaussian 98 style output file
      QFile g98(dir.filePath("g98.out"));
      QString g98Text;
      if (g98.open(QIODevice::ReadOnly)) g98Text = QString::fromLocal8Bit(g98.readAll());
      if (!parseG98(g98Text, input.atomicNumbers, result.vibrations)) {
         result.message = "Failed to read the vibrational frequencies from xtb\n" + tail(output);
         return result;
      }
      parseZeroPointEnergy(output, result.vibrations.zeroPointEnergy);
      result.hasVibrations = true;
      result.log += "--- g98.out ---\n" + g98Text;
   }

   if (input.orbitals && m_method != GFNFF) {
      QFile molden(dir.filePath("molden.input"));
      QString moldenText;
      if (molden.open(QIODevice::ReadOnly)) moldenText = QString::fromLocal8Bit(molden.readAll());
      if (parseMolden(moldenText, result.orbitals)) {
         result.hasOrbitals = true;
      }else {
         result.log += "--- the orbitals could not be read from molden.input ---\n";
      }
   }
   result.success = true;
   return result;
}

} } // end namespace IQmol::Calc
