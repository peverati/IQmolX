/*******************************************************************************

  This file is part of IQmol, a free molecular visualization program. See
  <http://iqmol.org> for more details.

  IQmol is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

********************************************************************************/

#include "XtbLibrary.h"
#include "XtbOutput.h"
#include "Optimizer.h"
#include "Vibrations.h"
#include <xtb.h>
#include <QCoreApplication>
#include <QElapsedTimer>
#include <QFile>
#include <QTemporaryDir>
#include <cmath>


namespace IQmol {
namespace Calc {

namespace {

// libxtb works in Bohr, IQmol in Angstrom
double const s_bohrPerAngstrom(1.0 / 0.529177210903);
double const s_debyePerEBohr(2.541746473);


// Owns the libxtb objects for one calculation
class Session {

   public:
      Session(Input const& input, XtbProcess::Method method) 
       : m_env(xtb_newEnvironment()), m_mol(0), m_calc(0), m_res(0), 
         m_natoms(input.atomicNumbers.size()), m_numbers(input.atomicNumbers.toStdVector())
      {
         xtb_setVerbosity(m_env, XTB_VERBOSITY_MUTED);

         std::vector<double> bohr(bohrCoordinates(input.coordinates));
         double charge(input.charge);
         int    unpaired(input.multiplicity - 1);

         m_mol  = xtb_newMolecule(m_env, &m_natoms, m_numbers.data(), bohr.data(), 
            &charge, &unpaired, 0, 0);
         m_calc = xtb_newCalculator();
         m_res  = xtb_newResults();
         if (failed()) return;

         switch (method) {
            case XtbProcess::GFN2:  xtb_loadGFN2xTB(m_env, m_mol, m_calc, 0); break;
            case XtbProcess::GFN1:  xtb_loadGFN1xTB(m_env, m_mol, m_calc, 0); break;
            case XtbProcess::GFNFF: xtb_loadGFNFF(m_env, m_mol, m_calc, 0);   break;
         }
         if (failed()) return;

         // GFN-FF has no electronic structure, so no electronic temperature either
         if (input.electronicTemperature >= 0.0 && method != XtbProcess::GFNFF) {
            xtb_setElectronicTemp(m_env, m_calc, input.electronicTemperature);
         }

         if (!input.solvent.isEmpty()) {
            // The solvent model has to be named, the library would otherwise use GBSA
            QByteArray solvent(input.solvent.toLatin1());
            char model[] = "alpb";
            xtb_setSolvent(m_env, m_calc, solvent.data(), 0, 0, 0, model);
         }
      }

      ~Session()
      {
         xtb_delResults(&m_res);
         xtb_delCalculator(&m_calc);
         xtb_delMolecule(&m_mol);
         xtb_delEnvironment(&m_env);
      }

      static std::vector<double> bohrCoordinates(QVector<double> const& angstrom)
      {
         std::vector<double> bohr(angstrom.size());
         for (int i = 0; i < angstrom.size(); ++i) bohr[i] = angstrom[i] * s_bohrPerAngstrom;
         return bohr;
      }

      bool failed() const { return xtb_checkEnvironment(m_env) != 0; }

      QString error() const
      {
         if (!failed()) return QString();
         char buffer[2048];
         int size(sizeof(buffer));
         xtb_getError(m_env, buffer, &size);
         return QString::fromLocal8Bit(buffer).trimmed();
      }

      // x in Bohr
      bool evaluate(std::vector<double> const& x, double& energy, std::vector<double>& gradient)
      {
         xtb_updateMolecule(m_env, m_mol, x.data(), 0);
         xtb_singlepoint(m_env, m_mol, m_calc, m_res);
         if (failed()) return false;
         xtb_getEnergy(m_env, m_res, &energy);
         gradient.assign(3 * m_natoms, 0.0);
         xtb_getGradient(m_env, m_res, gradient.data());
         return !failed();
      }

      // Numerical Hessian (Hartree/Bohr^2, 3N x 3N) and the derivatives of the dipole moment
      // with respect to the coordinates (in e, 3 values per coordinate).  The library wants a 
      // fresh results object here, as the Hessian is sensitive to the starting guess.
      bool hessian(std::vector<double> const& x, std::vector<double>& hessian, 
         std::vector<double>& dipoleGradient)
      {
         xtb_updateMolecule(m_env, m_mol, x.data(), 0);
         xtb_TResults fresh(xtb_newResults());
         hessian.assign(static_cast<size_t>(9) * m_natoms * m_natoms, 0.0);
         dipoleGradient.assign(static_cast<size_t>(9) * m_natoms, 0.0);
         // The xtb program converges the SCF more tightly for the displaced structures
         // of a Hessian (accuracy 0.3 instead of 1), which matters for soft modes
         xtb_setAccuracy(m_env, m_calc, 0.3);
         xtb_hessian(m_env, m_mol, m_calc, fresh, hessian.data(), 0, 0, dipoleGradient.data(), 0);
         xtb_setAccuracy(m_env, m_calc, 1.0);
         xtb_delResults(&fresh);
         return !failed();
      }

      // One more single point, this time recording xtb's printout and the dipole
      bool printout(std::vector<double> const& x, QString& text, double& dipoleDebye)
      {
         QTemporaryDir dir;
         if (!dir.isValid()) return false;
         QString path(dir.filePath("xtb.out"));

         xtb_setOutput(m_env, path.toLocal8Bit().constData());
         xtb_setVerbosity(m_env, XTB_VERBOSITY_FULL);
         double energy;
         std::vector<double> gradient;
         bool ok(evaluate(x, energy, gradient));
         xtb_releaseOutput(m_env);
         xtb_setVerbosity(m_env, XTB_VERBOSITY_MUTED);
         if (!ok) return false;

         QFile file(path);
         if (file.open(QIODevice::ReadOnly)) text = QString::fromLocal8Bit(file.readAll());

         double dipole[3] = { 0.0, 0.0, 0.0 };
         xtb_getDipole(m_env, m_res, dipole);
         dipoleDebye = s_debyePerEBohr * 
            std::sqrt(dipole[0]*dipole[0] + dipole[1]*dipole[1] + dipole[2]*dipole[2]);
         return true;
      }

   private:
      xtb_TEnvironment m_env;
      xtb_TMolecule    m_mol;
      xtb_TCalculator  m_calc;
      xtb_TResults     m_res;
      int              m_natoms;
      std::vector<int> m_numbers;
};


// libxtb does not check that the spin is consistent with the number of electrons
QString checkSpin(Input const& input)
{
   int electrons(-input.charge);
   for (int i = 0; i < input.atomicNumbers.size(); ++i) electrons += input.atomicNumbers[i];
   int unpaired(input.multiplicity - 1);

   if (electrons < 0 || unpaired > electrons || (electrons - unpaired) % 2 != 0) {
      return "Spin multiplicity " + QString::number(input.multiplicity) + " is inconsistent with " 
           + QString::number(electrons) + " electrons";
   }
   return QString();
}

} // end anonymous namespace


XtbLibrary::XtbLibrary(Method method) 
 : m_method(method), m_interactive(false), m_cancelled(false)
{
}


QString XtbLibrary::name() const
{
   switch (m_method) {
      case XtbProcess::GFN2:  return "GFN2-xTB";
      case XtbProcess::GFN1:  return "GFN1-xTB";
      case XtbProcess::GFNFF: return "GFN-FF";
   }
   return "xTB";
}


Result XtbLibrary::energy(Input const& input)
{
   return run(input, false);
}


Result XtbLibrary::optimize(Input const& input)
{
   return run(input, true);
}


Result XtbLibrary::frequencies(Input const& input)
{
   Result result;
   m_cancelled = false;

   int natoms(input.atomicNumbers.size());
   if (natoms == 0 || input.coordinates.size() != 3*natoms) {
      result.message = "Invalid input: inconsistent number of atoms and coordinates";
      return result;
   }
   // GFN-FF does not use the electrons, so any spin is acceptable
   QString spinError(m_method == XtbProcess::GFNFF ? QString() : checkSpin(input));
   if (!spinError.isEmpty()) {
      result.message = spinError;
      return result;
   }

   QElapsedTimer clock;
   clock.start();

   Session session(input, m_method);
   if (session.failed()) {
      result.message = "xtb setup failed\n" + session.error();
      return result;
   }

   std::vector<double> x(Session::bohrCoordinates(input.coordinates));
   QString log;
   auto report = [&](QString const& text) {
      log += text;
      if (m_progress) m_progress(text);
   };

   // The energy and gradient of the geometry: the frequencies only mean something at a 
   // stationary point
   QString text;
   double dipole(-1.0);
   if (!session.printout(x, text, dipole)) {
      result.message = "xtb failed\n" + session.error();
      return result;
   }
   report(text);
   if (!parseEnergy(text, result.energy)) {
      result.message = "Unable to read the energy from the xtb output";
      return result;
   }
   parseDetails(text, result);
   if (dipole >= 0.0) result.dipoleDebye = dipole;
   result.initialEnergy = result.energy;

   report(QString::asprintf("\nNumerical Hessian (%d displaced calculations) ...\n", 6 * natoms));
   if (m_interactive) QCoreApplication::processEvents();

   std::vector<double> hessian, dipoleGradient;
   if (!session.hessian(x, hessian, dipoleGradient)) {
      result.message = "xtb failed to calculate the Hessian\n" + session.error();
      result.log = log;
      return result;
   }

   std::vector<int> numbers(input.atomicNumbers.toStdVector());
   VibrationalAnalysis analysis(analyseVibrations(numbers, x, hessian, &dipoleGradient));
   if (!analysis.ok) {
      result.message = "Vibrational analysis failed: " + QString::fromStdString(analysis.message);
      result.log = log;
      return result;
   }

   VibrationalData& data(result.vibrations);
   for (size_t i = 0; i < analysis.frequencies.size(); ++i) {
       data.frequencies.append(analysis.frequencies[i]);
       data.reducedMasses.append(analysis.reducedMasses[i]);
       if (i < analysis.intensities.size()) data.intensities.append(analysis.intensities[i]);
       data.modes.append(QVector<double>::fromStdVector(analysis.modes[i]));
   }
   data.zeroPointEnergy = analysis.zeroPointEnergy;
   result.hasVibrations = true;

   report("\n  mode     frequency (cm-1)   IR intensity (km/mol)\n");
   for (int i = 0; i < data.frequencies.size(); ++i) {
       report(QString::asprintf("  %4d     %12.2f     %12.4f\n", i + 1, data.frequencies[i], 
          i < data.intensities.size() ? data.intensities[i] : 0.0));
   }
   report(QString::asprintf("\n  zero point energy  %.8f Eh\n", data.zeroPointEnergy));

   result.coordinates.resize(input.coordinates.size());
   for (int i = 0; i < input.coordinates.size(); ++i) result.coordinates[i] = input.coordinates[i];
   result.converged = true;
   result.seconds = clock.elapsed() / 1000.0;
   result.log = log;
   result.success = true;
   return result;
}


Result XtbLibrary::run(Input const& input, bool optimize)
{
   Result result;
   m_cancelled = false;

   int natoms(input.atomicNumbers.size());
   if (natoms == 0 || input.coordinates.size() != 3*natoms) {
      result.message = "Invalid input: inconsistent number of atoms and coordinates";
      return result;
   }
   if (input.multiplicity < 1) {
      result.message = "Invalid spin multiplicity";
      return result;
   }
   // GFN-FF does not use the electrons, so any spin is acceptable
   QString spinError(m_method == XtbProcess::GFNFF ? QString() : checkSpin(input));
   if (!spinError.isEmpty()) {
      result.message = spinError;
      return result;
   }

   QElapsedTimer clock;
   clock.start();

   Session session(input, m_method);
   if (session.failed()) {
      result.message = "xtb setup failed\n" + session.error();
      return result;
   }

   std::vector<double> x(Session::bohrCoordinates(input.coordinates));
   QString log;

   auto report = [&](QString const& text) {
      log += text;
      if (m_progress) m_progress(text);
   };

   if (optimize) {
      report("Geometry optimization with " + name() + " (libxtb) and IQmol's L-BFGS optimizer\n"
           "convergence: energy change < 5e-6 Eh, gradient norm < 1e-3 Eh/Bohr\n\n");

      bool first(true);
      auto evaluate = [&](std::vector<double> const& xx, double& e, std::vector<double>& g) {
         bool ok(session.evaluate(xx, e, g));
         if (ok && first) { result.initialEnergy = e; first = false; }
         return ok;
      };

      auto progress = [&](int cycle, double e, double gnorm, std::vector<double> const&) {
         report(QString::asprintf("  cycle %3d   E = %.8f Eh   |g| = %.6f Eh/Bohr\n", cycle, e, gnorm));
         if (m_interactive) QCoreApplication::processEvents();
         return !m_cancelled;
      };

      OptimizerResult opt(Optimizer().minimize(evaluate, x, progress));
      result.cycles = opt.cycles;
      result.converged = opt.converged;

      if (opt.cancelled) {
         result.cancelled = true;
         result.message = "Calculation cancelled";
         result.log = log;
         return result;
      }
      if (opt.x.empty() || (!opt.converged && session.failed())) {
         result.message = "xtb failed during the optimization\n" + session.error();
         result.log = log;
         return result;
      }

      x = opt.x;
      report(opt.converged ? "\n*** GEOMETRY OPTIMIZATION CONVERGED AFTER " 
                           + QString::number(opt.cycles) + " CYCLES ***\n\n"
                         : "\n*** " + QString::fromStdString(opt.message) + " ***\n\n");
      if (!opt.converged) result.message = "Geometry optimization did not converge\n" 
                                          + QString::fromStdString(opt.message);
   }else {
      result.converged = true;
   }

   // A final calculation at the result gives xtb's printout, the dipole and Gsolv
   QString text;
   double dipole(-1.0);
   if (!session.printout(x, text, dipole)) {
      result.message = "xtb failed\n" + session.error();
      result.log = log;
      return result;
   }
   report(text);

   if (!parseEnergy(text, result.energy)) {
      result.message = "Unable to read the energy from the xtb output";
      result.log = log;
      return result;
   }
   int cycles(result.cycles);
   double initialEnergy(result.initialEnergy);
   parseDetails(text, result);
   result.cycles = cycles;                       // the printout of a single point has no cycles
   if (optimize) result.initialEnergy = initialEnergy;
   if (dipole >= 0.0) result.dipoleDebye = dipole;

   result.coordinates.resize(x.size());
   for (size_t i = 0; i < x.size(); ++i) result.coordinates[i] = x[i] / s_bohrPerAngstrom;

   result.seconds = clock.elapsed() / 1000.0;
   result.log = log;
   result.success = true;
   return result;
}

} } // end namespace IQmol::Calc
