#pragma once
/*******************************************************************************

  This file is part of IQmol, a free molecular visualization program. See
  <http://iqmol.org> for more details.

  IQmol is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

********************************************************************************/

#include "Calculator.h"
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
         enum Method { GFN2, GFN1, GFNFF };

         explicit XtbProcess(QString const& executable = QString(), Method = GFN2);

         QString name() const override;

         void setExecutable(QString const& path) { m_executable = path; }
         void setMethod(Method method) { m_method = method; }
         void setTimeoutSeconds(int seconds) { m_timeout = seconds; }

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

         Result energy(Input const&) override;
         Result optimize(Input const&) override;
         Result frequencies(Input const&) override;

      private:
         enum Task { SinglePoint, Optimization, Hessian };
         Result run(Input const&, Task task);
         QStringList methodArguments() const;

         QString m_executable;
         Method  m_method;
         int     m_timeout;
         bool    m_interactive;
         bool    m_cancelled;
         QProcess* m_process;
         std::function<void(QString const&)> m_progress;
   };

} } // end namespace IQmol::Calc
