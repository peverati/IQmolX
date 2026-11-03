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

#include <QString>
#include <QStringList>
#include <functional>

class QProcess;

/*!
 *  \file PythonSetup.h
 *
 *  \brief Creates and extends the Python environment of the Python engines with uv
 *  (https://github.com/astral-sh/uv).  uv downloads a Python of its own, so no Python has
 *  to be installed on the computer; the environment lives in the application data folder 
 *  of the user.  Nothing runs without the user asking for it in the Python Engines window,
 *  and the commands that are run are written to the window as they start.
 *
 *  The steps are synchronous like the calculations: with setInteractive(true) the Qt event
 *  loop keeps running meanwhile, so the window stays responsive and cancel() can be called.
 */

namespace IQmol {
namespace Calc {

   class PythonSetup {
      public:
         PythonSetup();

         /// Where uv is looked for: the hint, $IQMOLX_UV, the Resources folder of the 
         /// application, the PATH and the usual install folders.  Empty if it was not found.
         static QString findUv(QString const& hint = QString());

         /// The folder of the environment that IQmol makes: in the application data folder
         static QString defaultEnvironment();

         /// The interpreter of an environment folder (empty if there is none)
         static QString pythonIn(QString const& environment);

         /// The folder of the iqmolx-engines package (with its pyproject.toml): from
         /// $IQMOLX_ENGINES_PATH (its src folder) or the Resources folder of the application.
         static QString helperSource();

         void setUv(QString const& path) { m_uv = path; }
         void setInteractive(bool interactive) { m_interactive = interactive; }
         void setProgressCallback(std::function<void(QString const&)> callback) 
         { 
            m_progress = callback; 
         }
         void cancel() { m_cancelled = true; }

         /// Makes the environment (a Python that uv downloads if necessary) and installs the
         /// helper into it.  Returns false and the reason in 'error' on failure.
         bool createEnvironment(QString const& environment, QString& error);

         /// Installs pip requirements (such as "torch") into the environment of an interpreter
         bool install(QString const& python, QStringList const& requirements, QString& error);

         bool cancelled() const { return m_cancelled; }

      private:
         bool run(QStringList const& arguments, QString& error);

         QString m_uv;
         bool    m_interactive;
         bool    m_cancelled;
         std::function<void(QString const&)> m_progress;
   };

} } // end namespace IQmol::Calc
