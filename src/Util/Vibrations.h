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

#include <string>
#include <vector>

/*!
 *  \file Vibrations.h
 *
 *  \brief Harmonic vibrational analysis: turns a Cartesian Hessian (and optionally the
 *  derivatives of the dipole moment) into frequencies, normal modes and IR intensities.
 *  It has no dependencies on Qt or on the program that produced the Hessian.
 */

namespace IQmol {
namespace Calc {

   struct VibrationalAnalysis {
      bool   ok;
      std::string message;

      /// Frequencies in cm^-1, ascending.  Imaginary frequencies are reported as negative.
      std::vector<double> frequencies;

      /// IR intensities in km/mol, empty if no dipole derivatives were given
      std::vector<double> intensities;

      /// Reduced masses in amu
      std::vector<double> reducedMasses;

      /// Normal modes: for each mode the Cartesian displacement of every atom
      /// (3N numbers, x y z per atom), normalized so that the squares sum to one
      std::vector<std::vector<double> > modes;

      /// Zero point energy in Hartree, from the real frequencies
      double zeroPointEnergy;

      /// Number of translations and rotations that were removed (6, or 5 for a linear molecule)
      int translationsRotations;

      VibrationalAnalysis() : ok(false), zeroPointEnergy(0.0), translationsRotations(0) { }
   };

   /// Atomic mass in amu (the same table that xtb uses), 0 if the element is unknown
   double atomicMass(int atomicNumber);

   /// Harmonic analysis.
   ///   atomicNumbers   N atoms
   ///   bohr            3N coordinates in Bohr
   ///   hessian         3N x 3N second derivatives of the energy in Hartree/Bohr^2 (row major)
   ///   dipoleGradient  optional, derivative of each dipole component with respect to each
   ///                   coordinate in e: for coordinate j the 3 values dmu_x/dx_j, dmu_y/dx_j,
   ///                   dmu_z/dx_j follow each other, 9N numbers in all
   VibrationalAnalysis analyseVibrations(std::vector<int> const& atomicNumbers,
      std::vector<double> const& bohr, std::vector<double> const& hessian,
      std::vector<double> const* dipoleGradient = 0);

   /// Eigenvalues (ascending) and eigenvectors (the columns of the row major n x n matrix 
   /// 'vectors') of a real symmetric n x n matrix, by Jacobi rotations.  The matrix is destroyed.
   void symmetricEigenvalues(std::vector<double>& matrix, int n, std::vector<double>& values,
      std::vector<double>& vectors);

} } // end namespace IQmol::Calc
