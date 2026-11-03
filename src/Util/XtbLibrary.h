#pragma once
/*******************************************************************************

  Copyright (C) 2026 Roberto Peverati

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

/*!
 *  \file XtbLibrary.h
 *
 *  \brief Calculator that runs xtb inside IQmol by calling libxtb, instead of
 *  starting the xtb program.  libxtb provides energies and gradients but no
 *  geometry optimizer, so optimizations use IQmol's own L-BFGS optimizer.
 *
 *  Only compiled when IQmol is built with XTB_ROOT (IQMOL_WITH_LIBXTB).
 */

namespace IQmol {
namespace Calc {

   class XtbLibrary : public Calculator {
      public:
         typedef XtbProcess::Method Method;

         explicit XtbLibrary(Method = XtbProcess::GFN2);

         QString name() const override;

         void setMethod(Method method) { m_method = method; }

         /// libxtb does not provide the basis the orbitals are expanded in, so when orbitals
         /// are requested the xtb program is run once more at the final geometry to get them.
         void setExecutable(QString const& path) { m_executable = path; }
         void setProgressCallback(std::function<void(QString const&)> callback) override 
         { 
            m_progress = callback; 
         }
         void setInteractive(bool interactive) override { m_interactive = interactive; }
         void cancel() override;

         Result energy(Input const&) override;
         Result optimize(Input const&) override;
         Result frequencies(Input const&) override;

      private:
         Result run(Input const&, bool optimize);
         Result vibrations(Input const&);
         Result withOrbitals(Input const&, Result);

         QString m_executable;
         XtbProcess* m_process;
         Method  m_method;
         bool    m_interactive;
         bool    m_cancelled;
         std::function<void(QString const&)> m_progress;
   };

} } // end namespace IQmol::Calc
