#pragma once

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

#include "Calculator.h"
#include <QJsonArray>
#include <QJsonObject>
#include <QStringList>
#include <functional>
#include <vector>

/*!
 *  \file PythonEngine.h
 *
 *  \brief Calculator that uses an engine of the iqmolx-engines Python helper (machine
 *  learning potentials, ASE calculators, ...).  The helper is a child process that stays
 *  alive between calculations, so a model is loaded once.  IQmol talks to it with one JSON
 *  message per line, as described in python/iqmolx-engines/PROTOCOL.md; nothing here knows
 *  about a specific engine, everything comes from the capabilities the helper reports.
 *
 *  The helper works in Angstrom and eV, IQmol in Bohr and Hartree.  The helper provides
 *  energies and forces (and perhaps a Hessian); geometry optimizations use IQmol's L-BFGS
 *  optimizer and frequencies use Vibrations, like the libxtb backend.
 */

namespace IQmol {
namespace Calc {

   class PythonHelper;
   struct PythonFailure;

   class PythonEngine : public Calculator {
      public:
         explicit PythonEngine(QString const& python = QString(), 
            QString const& engine = QString());
         ~PythonEngine();

         /// The title of the engine once its capabilities are known, else its name
         QString name() const override;

         /// The Python interpreter that runs the helper (default: python3 from the PATH)
         void setPython(QString const& path);

         /// Directory added to PYTHONPATH, for a helper that is not installed in the
         /// environment of the interpreter (the source folder of iqmolx-engines)
         void setPythonPath(QString const& path);

         /// The engine to use (a name from listEngines) and its settings.  The load options
         /// are applied when the model loads, so changing them loads the model again; the 
         /// compute options go with every calculation.  Both are objects of option name to value.
         void setEngine(QString const& name);
         void setLoadOptions(QJsonObject const& options);
         void setComputeOptions(QJsonObject const& options);

         /// Use the Hessian of the engine when it offers one (the default).  Otherwise, or 
         /// when it does not, the Hessian is built from forces at displaced geometries.
         void setAnalyticHessian(bool use) { m_analyticHessian = use; }

         /// Displacement of the numerical Hessian in Bohr
         void setHessianStep(double bohr) { m_hessianStep = bohr; }

         /// The engines of the helper with their capabilities (the engine descriptors of the 
         /// protocol).  Starts the helper if necessary; empty on failure, with the reason in 
         /// 'error'.
         QJsonArray listEngines(QString* error = 0);

         /// The version of the helper package that is running, empty if it is not running
         QString helperVersion() const;

         /// The descriptor of the selected engine, empty if it is not (yet) known
         QJsonObject engineInfo() const { return m_descriptor; }

         /// Stops the helper process.  It starts again with the next calculation.
         void shutdown();

         void setProgressCallback(std::function<void(QString const&)> callback) override 
         { 
            m_progress = callback; 
         }
         void setInteractive(bool interactive) override { m_interactive = interactive; }
         void cancel() override { m_cancelled = true; }

         Result energy(Input const&) override;
         Result optimize(Input const&) override;
         Result frequencies(Input const&) override;

      private:
         bool prepare(PythonFailure&);
         bool compute(Input const&, std::vector<double> const& bohr, QStringList const& want,
            QJsonObject& result, PythonFailure&);
         bool gradient(Input const&, std::vector<double> const& bohr, double& energy, 
            std::vector<double>& gradient, PythonFailure&, QJsonObject* details = 0);
         bool hessian(Input const&, std::vector<double> const& bohr, 
            std::vector<double>& hessian, PythonFailure&);
         bool offers(QString const& property) const;
         void report(QString const& text);
         Result failed(PythonFailure const&, QString const& log);
         static bool validInput(Input const&, Result&);

         QString m_python;
         QString m_pythonPath;
         QString m_engine;
         QJsonObject m_loadOptions;
         QJsonObject m_computeOptions;
         QJsonObject m_descriptor;
         QJsonArray  m_engines;
         QString m_loadedKey;           // engine and load options of the model in the helper
         PythonHelper* m_helper;
         bool   m_analyticHessian;
         double m_hessianStep;
         bool   m_interactive;
         bool   m_cancelled;
         QString m_log;
         std::function<void(QString const&)> m_progress;
   };

} } // end namespace IQmol::Calc
