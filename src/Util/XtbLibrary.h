#pragma once
/*******************************************************************************

  This file is part of IQmol, a free molecular visualization program. See
  <http://iqmol.org> for more details.

  IQmol is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

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
         void setProgressCallback(std::function<void(QString const&)> callback) override 
         { 
            m_progress = callback; 
         }
         void setInteractive(bool interactive) override { m_interactive = interactive; }
         void cancel() override { m_cancelled = true; }

         Result energy(Input const&) override;
         Result optimize(Input const&) override;

      private:
         Result run(Input const&, bool optimize);

         Method  m_method;
         bool    m_interactive;
         bool    m_cancelled;
         std::function<void(QString const&)> m_progress;
   };

} } // end namespace IQmol::Calc
