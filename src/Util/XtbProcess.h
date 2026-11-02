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

         QString name() const;

         void setExecutable(QString const& path) { m_executable = path; }
         void setMethod(Method method) { m_method = method; }
         void setTimeoutSeconds(int seconds) { m_timeout = seconds; }

         /// Looks for xtb: the given hint, then $XTB_PATH, then PATH.
         /// Returns an empty string if it cannot be found.
         static QString findExecutable(QString const& hint = QString());

         Result energy(Input const&);
         Result optimize(Input const&);

      private:
         Result run(Input const&, bool optimize);
         QStringList methodArguments() const;

         QString m_executable;
         Method  m_method;
         int     m_timeout;
   };

} } // end namespace IQmol::Calc
