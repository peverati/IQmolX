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

   /// Zero point energy of a frequency calculation, in Hartree
   bool parseZeroPointEnergy(QString const& output, double& energy);

   /// Reads the frequencies, IR intensities, reduced masses and normal modes from the 
   /// g98.out file that xtb writes for a Hessian calculation (Gaussian 98 format) for a
   /// molecule with these atoms.  xtb writes the normal modes mass weighted; they are 
   /// converted to Cartesian displacements and normalized.
   bool parseG98(QString const& text, QVector<int> const& atomicNumbers, VibrationalData& data);

   /// Reads the molecular orbitals from the molden.input file that xtb writes with --molden.
   /// xtb gives the contraction coefficients for unnormalized primitives, they are converted
   /// to those of normalized primitives here.  Only the first set of orbitals (alpha) is read;
   /// xtb uses the same spatial orbitals for both spins.  With tbliteConvention the file is
   /// that of g-xTB, whose Cartesian d functions are normalized differently.
   bool parseMolden(QString const& text, OrbitalData& data, bool tbliteConvention = false);

} } // end namespace IQmol::Calc
