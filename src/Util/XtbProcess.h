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
#include <QStringList>
#include <functional>

class QProcess;

/*!
 *  \file XtbProcess.h
 *
 *  \brief Calculator that runs the external xtb program (Grimme group) in a
 *  temporary directory.  A library-based implementation can later replace this
 *  behind the same Calculator interface.
 */

namespace IQmol {
namespace Calc {

   class XtbProcess : public Calculator {
      public:
         enum Method { GFN2, GFN1, GFNFF, GXTB };

         explicit XtbProcess(QString const& executable = QString(), Method = GFN2);

         QString name() const override;

         void setExecutable(QString const& path) { m_executable = path; }
         void setMethod(Method method) { m_method = method; }

         /// The implicit solvation model of g-xTB, "gbe" (generalized Born) or "cosmo" (ddCOSMO).
         /// The ALPB model of the other methods is not available for g-xTB.
         void setSolvationModel(QString const& model) { m_solvationModel = model; }
         void setTimeoutSeconds(int seconds) { m_timeout = seconds; }

         /// Further command line options, added after those IQmol sets so they can override them
         void setExtraArguments(QStringList const& arguments) { m_extraArguments = arguments; }

         /// Directory to run in.  The input and output files stay there (an existing
         /// xtbopt.xyz, molden.input etc. is removed first).  By default a temporary
         /// directory is used and removed afterwards.
         void setWorkingDirectory(QString const& path) { m_workingDirectory = path; }

         /// Called with each chunk of xtb output as it arrives.  This only has an effect
         /// together with setInteractive(true).
         void setProgressCallback(std::function<void(QString const&)> callback) override
         { 
            m_progress = callback; 
         }

         /// When true the Qt event loop keeps running while xtb does, so a GUI stays
         /// responsive (and cancel() can be called).  The default is to block.
         void setInteractive(bool interactive) override { m_interactive = interactive; }

         /// Kills a running calculation (interactive mode only)
         void cancel() override;

         /// Looks for xtb: the given hint, then $XTB_PATH, then PATH.
         /// Returns an empty string if it cannot be found.
         static QString findExecutable(QString const& hint = QString());

         enum Task { SinglePoint, Optimization, Hessian };

         /// By default the temporary directory is removed after each calculation.  When this
         /// is set it is kept (and reused by the following calculations of this object) and
         /// reported in Result::directory, whose owner then has to remove it.
         void setKeepTemporaryDirectory(bool keep) { m_keepTemporary = keep; }

         Result energy(Input const&) override;
         Result optimize(Input const&) override;
         Result frequencies(Input const&) override;

      private:
         Result run(Input const&, Task task);
         QStringList arguments(Input const&, Task) const;
         QStringList methodArguments() const;

         QString m_executable;
         QStringList m_extraArguments;
         QString m_workingDirectory;
         bool m_keepTemporary;
         QString m_temporaryPath;
         Method  m_method;
         QString m_solvationModel;
         int     m_timeout;
         bool    m_interactive;
         bool    m_cancelled;
         QProcess* m_process;
         std::function<void(QString const&)> m_progress;
   };

} } // end namespace IQmol::Calc
