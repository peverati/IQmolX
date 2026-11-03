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
 *  \file XtbOutput.h
 *
 *  \brief Reads values from the text output of xtb.  Shared by the backend that
 *  runs the xtb program and the one that calls libxtb, so both report the same
 *  quantities in the same way.
 */

namespace IQmol {
namespace Calc {

   /// Last 'TOTAL ENERGY' in the output, in Hartree
   bool parseEnergy(QString const& output, double& energy);

   /// Fills in the optional details of the result (gradient norm, HOMO-LUMO gap, 
   /// number of cycles, dipole, solvation free energy, energy of the first geometry)
   void parseDetails(QString const& output, Result& result);

} } // end namespace IQmol::Calc
