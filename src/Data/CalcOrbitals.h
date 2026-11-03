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

#include "Data/CanonicalOrbitals.h"
#include "Util/Calculator.h"
#include "QGLViewer/vec.h"

namespace IQmol {
namespace Data {

   /// Converts the orbitals of a Calc::Result (bohr exponents, one set of spatial orbitals)
   /// to IQmol's canonical orbitals.  The positions of the atoms are in Angstrom.  The
   /// numbers of alpha and beta electrons follow from the occupations (the valence electrons
   /// of the xTB methods) and the spin multiplicity.  Returns 0 if the data are inconsistent.
   CanonicalOrbitals* makeCanonicalOrbitals(Calc::OrbitalData const&, 
      QList<qglviewer::Vec> const& positions, int multiplicity, QString const& label);

} } // end namespace IQmol::Data
