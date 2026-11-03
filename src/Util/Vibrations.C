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

#include "Vibrations.h"
#include <algorithm>
#include <cmath>


namespace IQmol {
namespace Calc {

namespace {

// sqrt(Eh / (a0^2 amu)) / (2 pi c) in cm^-1
double const s_frequencyFactor(5140.487446);

// Converts (dmu/dQ)^2 in e^2/amu to the IR intensity in km/mol:  N_A pi / (3 c^2 4 pi eps0)
double const s_intensityFactor(974.8801);

// 1 cm^-1 in Hartree
double const s_wavenumberToHartree(4.556335252912e-6);

// Atomic masses in amu, as in xtb (mctc-lib, standard atomic weights)
double const s_masses[118] = {
      1.00794075, 4.00260193, 6.94003660, 9.01218307, 10.81102805, 12.01073590,
      14.00670321, 15.99940492, 18.99840316, 20.18004638, 22.98976928, 24.30505162,
      26.98153853, 28.08549871, 30.97376200, 32.06478741, 35.45293758, 39.94779856,
      39.09830091, 40.07802251, 44.95590828, 47.86674496, 50.94146504, 51.99613176,
      54.93804391, 55.84514443, 58.93319429, 58.69334711, 63.54603995, 65.37778253,
      69.72306607, 72.62755016, 74.92159457, 78.95938856, 79.90352778, 83.79800000,
      85.46766360, 87.61664447, 88.90584030, 91.22364160, 92.90637300, 95.95978854,
      97.90721240, 101.06494014, 102.90549800, 106.41532751, 107.86814963, 112.41155782,
      114.81808663, 118.71011259, 121.75978367, 127.60312648, 126.90447190, 131.29276145,
      132.90545196, 137.32689163, 138.90546887, 140.11573074, 140.90765760, 144.24159603,
      144.91275590, 150.36635571, 151.96437813, 157.25213065, 158.92535470, 162.49947282,
      164.93032880, 167.25908265, 168.93421790, 173.05415017, 174.96681496, 178.48497872,
      180.94787564, 183.84177755, 186.20670455, 190.22485963, 192.21605165, 195.08445686,
      196.96656879, 200.59916703, 204.38341284, 207.21690806, 208.98039910, 208.98243080,
      209.98714790, 222.01757820, 223.01973600, 226.02541030, 227.02775230, 232.03805580,
      231.03588420, 238.02891046, 237.04817360, 244.06420530, 243.06138130, 247.07035410,
      247.07030730, 251.07958860, 252.08298000, 257.09510610, 258.09843150, 259.10103000,
      262.10961000, 267.12179000, 269.12791000, 271.13393000, 270.13336000, 276.14846000,
      276.15159000, 280.16131000, 282.16912000, 284.17416000, 284.17873000, 289.19042000,
      288.19274000, 293.20449000, 292.20746000, 294.21392000
};


double dot(std::vector<double> const& a, std::vector<double> const& b)
{
   double sum(0.0);
   for (size_t i = 0; i < a.size(); ++i) sum += a[i] * b[i];
   return sum;
}

} // end anonymous namespace


double atomicMass(int atomicNumber)
{
   return (atomicNumber >= 1 && atomicNumber <= 118) ? s_masses[atomicNumber - 1] : 0.0;
}


void symmetricEigenvalues(std::vector<double>& a, int n, std::vector<double>& values, 
   std::vector<double>& vectors)
{
   vectors.assign(static_cast<size_t>(n) * n, 0.0);
   for (int i = 0; i < n; ++i) vectors[i * n + i] = 1.0;

   for (int sweep = 0; sweep < 100; ++sweep) {
       double off(0.0), diagonal(0.0);
       for (int i = 0; i < n; ++i) {
           diagonal += a[i * n + i] * a[i * n + i];
           for (int j = i + 1; j < n; ++j) off += a[i * n + j] * a[i * n + j];
       }
       if (off <= 1e-26 * (diagonal + 1e-300) || off < 1e-300) break;

       for (int p = 0; p < n - 1; ++p) {
           for (int q = p + 1; q < n; ++q) {
               double apq(a[p * n + q]);
               if (std::fabs(apq) < 1e-300) continue;

               double theta((a[q * n + q] - a[p * n + p]) / (2.0 * apq));
               double t((theta >= 0.0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0)));
               double c(1.0 / std::sqrt(t * t + 1.0));
               double s(t * c);

               for (int k = 0; k < n; ++k) {         // columns p and q
                   double akp(a[k * n + p]), akq(a[k * n + q]);
                   a[k * n + p] = c * akp - s * akq;
                   a[k * n + q] = s * akp + c * akq;
               }
               for (int k = 0; k < n; ++k) {         // rows p and q
                   double apk(a[p * n + k]), aqk(a[q * n + k]);
                   a[p * n + k] = c * apk - s * aqk;
                   a[q * n + k] = s * apk + c * aqk;
               }
               for (int k = 0; k < n; ++k) {         // eigenvectors
                   double vkp(vectors[k * n + p]), vkq(vectors[k * n + q]);
                   vectors[k * n + p] = c * vkp - s * vkq;
                   vectors[k * n + q] = s * vkp + c * vkq;
               }
           }
       }
   }

   // Sort ascending
   std::vector<int> order(n);
   for (int i = 0; i < n; ++i) order[i] = i;
   std::sort(order.begin(), order.end(), [&a, n](int x, int y) { return a[x * n + x] < a[y * n + y]; });

   values.resize(n);
   std::vector<double> sorted(static_cast<size_t>(n) * n);
   for (int k = 0; k < n; ++k) {
       values[k] = a[order[k] * n + order[k]];
       for (int i = 0; i < n; ++i) sorted[i * n + k] = vectors[i * n + order[k]];
   }
   vectors.swap(sorted);
}


VibrationalAnalysis analyseVibrations(std::vector<int> const& atomicNumbers,
   std::vector<double> const& bohr, std::vector<double> const& hessian,
   std::vector<double> const* dipoleGradient)
{
   VibrationalAnalysis result;
   int const natoms(static_cast<int>(atomicNumbers.size()));
   int const n(3 * natoms);

   if (natoms == 0 || bohr.size() != static_cast<size_t>(n) || 
       hessian.size() != static_cast<size_t>(n) * n) {
      result.message = "Inconsistent sizes of the coordinates and the Hessian";
      return result;
   }
   if (dipoleGradient && dipoleGradient->size() != static_cast<size_t>(3) * n) {
      result.message = "Inconsistent size of the dipole derivatives";
      return result;
   }

   std::vector<double> mass(natoms), sqrtMass(natoms);
   for (int i = 0; i < natoms; ++i) {
       mass[i] = atomicMass(atomicNumbers[i]);
       if (mass[i] <= 0.0) {
          result.message = "Unknown element, atomic number " + std::to_string(atomicNumbers[i]);
          return result;
       }
       sqrtMass[i] = std::sqrt(mass[i]);
   }

   // Mass weighted Hessian, symmetrized
   std::vector<double> weighted(static_cast<size_t>(n) * n);
   for (int i = 0; i < n; ++i) {
       for (int j = 0; j < n; ++j) {
           double h(0.5 * (hessian[i * n + j] + hessian[j * n + i]));
           weighted[i * n + j] = h / (sqrtMass[i / 3] * sqrtMass[j / 3]);
       }
   }

   // Translations and rotations about the centre of mass, orthonormalized.  A linear 
   // molecule has only two rotations and a single atom none, which shows as a vector
   // with no length left.
   double totalMass(0.0), centre[3] = { 0.0, 0.0, 0.0 };
   for (int i = 0; i < natoms; ++i) {
       totalMass += mass[i];
       for (int k = 0; k < 3; ++k) centre[k] += mass[i] * bohr[3 * i + k];
   }
   for (int k = 0; k < 3; ++k) centre[k] /= totalMass;

   std::vector<std::vector<double> > external;
   auto addExternal = [&external, n](std::vector<double> v) {
      double reference(std::sqrt(dot(v, v)));
      if (reference < 1e-12) return;
      for (int pass = 0; pass < 2; ++pass) {
          for (size_t t = 0; t < external.size(); ++t) {
              double overlap(dot(v, external[t]));
              for (int i = 0; i < n; ++i) v[i] -= overlap * external[t][i];
          }
      }
      double length(std::sqrt(dot(v, v)));
      if (length < 1e-6 * reference || length < 1e-8) return;
      for (int i = 0; i < n; ++i) v[i] /= length;
      external.push_back(v);
   };

   for (int axis = 0; axis < 3; ++axis) {
       std::vector<double> v(n, 0.0);
       for (int i = 0; i < natoms; ++i) v[3 * i + axis] = sqrtMass[i];
       addExternal(v);
   }
   for (int axis = 0; axis < 3; ++axis) {
       std::vector<double> v(n, 0.0);
       int a((axis + 1) % 3), b((axis + 2) % 3);
       for (int i = 0; i < natoms; ++i) {
           double ra(bohr[3 * i + a] - centre[a]), rb(bohr[3 * i + b] - centre[b]);
           // e_axis x r = (0,..,-rb, ra,..) in the components a and b
           v[3 * i + a] = -sqrtMass[i] * rb;
           v[3 * i + b] =  sqrtMass[i] * ra;
       }
       addExternal(v);
   }
   result.translationsRotations = static_cast<int>(external.size());

   int const m(n - result.translationsRotations);
   if (m <= 0) { result.ok = true; return result; }   // a single atom has no vibrations

   // Orthonormal basis of the space orthogonal to the translations and rotations, built from
   // the unit vectors with the largest remainder first (Gram-Schmidt with pivoting)
   std::vector<std::vector<double> > remainder(n, std::vector<double>(n, 0.0));
   for (int i = 0; i < n; ++i) {
       remainder[i][i] = 1.0;
       for (size_t t = 0; t < external.size(); ++t) {
           double overlap(external[t][i]);
           for (int k = 0; k < n; ++k) remainder[i][k] -= overlap * external[t][k];
       }
   }
   std::vector<std::vector<double> > basis;
   for (int step = 0; step < m; ++step) {
       int best(-1);
       double bestLength(-1.0);
       for (int i = 0; i < n; ++i) {
           double length(dot(remainder[i], remainder[i]));
           if (length > bestLength) { bestLength = length; best = i; }
       }
       std::vector<double> b(remainder[best]);
       double length(std::sqrt(bestLength));
       if (length < 1e-8) { result.message = "Unable to construct the vibrational space"; return result; }
       for (int k = 0; k < n; ++k) b[k] /= length;
       for (int i = 0; i < n; ++i) {
           double overlap(dot(remainder[i], b));
           for (int k = 0; k < n; ++k) remainder[i][k] -= overlap * b[k];
       }
       basis.push_back(b);
   }

   // Hessian in that space, and its eigenvalues (Hartree / (Bohr^2 amu))
   std::vector<double> reduced(static_cast<size_t>(m) * m, 0.0);
   {
      std::vector<std::vector<double> > hb(m, std::vector<double>(n, 0.0));
      for (int a = 0; a < m; ++a) {
          for (int i = 0; i < n; ++i) {
              double sum(0.0);
              for (int j = 0; j < n; ++j) sum += weighted[i * n + j] * basis[a][j];
              hb[a][i] = sum;
          }
      }
      for (int a = 0; a < m; ++a) {
          for (int b = a; b < m; ++b) {
              double sum(dot(basis[a], hb[b]));
              reduced[a * m + b] = sum;
              reduced[b * m + a] = sum;
          }
      }
   }

   std::vector<double> values, vectors;
   symmetricEigenvalues(reduced, m, values, vectors);

   bool const haveIntensities(dipoleGradient != 0);
   for (int k = 0; k < m; ++k) {
       double lambda(values[k]);
       double frequency(s_frequencyFactor * std::sqrt(std::fabs(lambda)));
       if (lambda < 0.0) frequency = -frequency;

       // Mass weighted normal mode in the full space, then Cartesian displacements
       std::vector<double> v(n, 0.0);
       for (int a = 0; a < m; ++a) {
           double w(vectors[a * m + k]);
           for (int i = 0; i < n; ++i) v[i] += w * basis[a][i];
       }
       std::vector<double> d(n);
       double sumD(0.0), inverseMass(0.0);
       for (int i = 0; i < n; ++i) {
           d[i] = v[i] / sqrtMass[i / 3];
           sumD += d[i] * d[i];
           inverseMass += v[i] * v[i] / mass[i / 3];
       }

       if (haveIntensities) {
          double derivative2(0.0);
          for (int comp = 0; comp < 3; ++comp) {
              double sum(0.0);
              for (int j = 0; j < n; ++j) sum += (*dipoleGradient)[3 * j + comp] * d[j];
              derivative2 += sum * sum;
          }
          result.intensities.push_back(s_intensityFactor * derivative2);
       }

       double norm(std::sqrt(sumD));
       for (int i = 0; i < n; ++i) d[i] /= norm;

       result.frequencies.push_back(frequency);
       result.reducedMasses.push_back(1.0 / inverseMass);
       result.modes.push_back(d);
       if (frequency > 0.0) result.zeroPointEnergy += 0.5 * frequency * s_wavenumberToHartree;
   }

   result.ok = true;
   return result;
}

} } // end namespace IQmol::Calc
