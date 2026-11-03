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

#include "CalcOrbitals.h"
#include "Util/Constants.h"
#include <cmath>

namespace IQmol {
namespace Data {

CanonicalOrbitals* makeCanonicalOrbitals(Calc::OrbitalData const& orbitals, 
   QList<qglviewer::Vec> const& positions, int multiplicity, QString const& label)
{
   if (orbitals.coefficients.isEmpty() || orbitals.energies.size() != orbitals.coefficients.size()) {
      return 0;
   }

   // The exponents are in bohr^-2, IQmol works in Angstrom
   static double const toAngstrom(std::pow(Constants::BohrToAngstrom, -2.0));
   static Shell::AngularMomentum const momentum[3] = { Shell::S, Shell::P, Shell::D6 };

   ShellList shells;
   for (Calc::OrbitalShell const& shell : orbitals.shells) {
       if (shell.atom < 0 || shell.atom >= positions.size() || shell.angularMomentum < 0
          || shell.angularMomentum > 2) {
          return 0;
       }
       QList<double> exponents, coefficients;
       for (int i = 0; i < shell.exponents.size(); ++i) {
           exponents.append(shell.exponents[i] * toAngstrom);
           coefficients.append(shell.coefficients[i]);
       }
       shells.append(new Shell(momentum[shell.angularMomentum], shell.atom, 
          positions[shell.atom], exponents, coefficients));
   }

   // The occupations add up to the number of electrons also with a smeared occupation
   double occupation(0.0);
   for (double n : orbitals.occupations) occupation += n;
   int electrons(qRound(occupation));
   int nBeta((electrons - (multiplicity - 1)) / 2);
   int nAlpha(electrons - nBeta);
   if (nBeta < 0 || nAlpha > orbitals.coefficients.size()) return 0;

   QList<double> coefficients, energies;
   for (int i = 0; i < orbitals.coefficients.size(); ++i) {
       coefficients << orbitals.coefficients[i].toList();
       energies.append(orbitals.energies[i]);
   }

   CanonicalOrbitals* canonical(new CanonicalOrbitals(nAlpha, nBeta, shells, coefficients, 
      energies, QList<double>(), QList<double>(), label));
   if (!canonical->consistent()) {
      delete canonical;
      return 0;
   }
   return canonical;
}

} }
