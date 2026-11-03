// Orbitals from the xtb program: reads molden.input and checks numerically, on a grid,
// that the orbitals are orthonormal in IQmol's basis function conventions (normalization of
// the primitives, ordering and scaling of the p and Cartesian d functions).
// Usage: XtbOrbitalsTest path/to/xtb

#include "Data/CalcOrbitals.h"
#include "Util/XtbProcess.h"
#include "Math/Matrix.h"
#include <QCoreApplication>
#include <cmath>
#include <cstdio>

using namespace IQmol;

static int s_failures(0);

static void check(bool ok, char const* what)
{
   std::printf("%s: %s\n", ok ? "PASS" : "FAIL", what);
   if (!ok) ++s_failures;
}

// returns the largest deviation of <i|j> from the unit matrix for the lowest nCheck orbitals
static double orthonormality(Calc::Result const& r, int nCheck, unsigned nAlpha, int multiplicity = 1)
{
   QList<qglviewer::Vec> positions;
   for (int i = 0; i < r.coordinates.size() / 3; ++i) {
       positions.append(qglviewer::Vec(r.coordinates[3*i], r.coordinates[3*i+1], r.coordinates[3*i+2]));
   }
   Data::CanonicalOrbitals* orbs(Data::makeCanonicalOrbitals(r.orbitals, positions, multiplicity, "test"));
   if (!orbs) return 1e9;

   Data::ShellList& shells(orbs->shellList());
   Matrix const& C(orbs->alphaCoefficients());
   int nBasis(orbs->nBasis());

   double lo[3] = {1e9, 1e9, 1e9}, hi[3] = {-1e9, -1e9, -1e9};
   for (auto const& p : positions) {
       for (int c = 0; c < 3; ++c) { lo[c] = std::min(lo[c], p[c]); hi[c] = std::max(hi[c], p[c]); }
   }
   double const margin(6.0), h(0.045);
   std::vector<double> ao(nBasis), values;
   std::vector<double> S(nCheck*nCheck, 0.0), mo(nCheck);
   int n[3];
   for (int c = 0; c < 3; ++c) n[c] = (int)((hi[c]-lo[c]+2*margin)/h) + 1;

   for (int ix = 0; ix < n[0]; ++ix) {
      double x(lo[0]-margin+ix*h);
      for (int iy = 0; iy < n[1]; ++iy) {
         double y(lo[1]-margin+iy*h);
         for (int iz = 0; iz < n[2]; ++iz) {
            double z(lo[2]-margin+iz*h);
            int k(0);
            for (int s = 0; s < shells.size(); ++s) {
                unsigned nf(shells[s]->nBasis());
                if (shells[s]->evaluate(x, y, z, values)) {
                   for (unsigned f = 0; f < nf; ++f) ao[k+f] = values[f];
                }else {
                   for (unsigned f = 0; f < nf; ++f) ao[k+f] = 0.0;
                }
                k += nf;
            }
            for (int i = 0; i < nCheck; ++i) {
                double sum(0.0);
                for (int j = 0; j < nBasis; ++j) sum += C(i,j)*ao[j];
                mo[i] = sum;
            }
            for (int i = 0; i < nCheck; ++i) for (int j = 0; j <= i; ++j) S[i*nCheck+j] += mo[i]*mo[j];
         }
      }
   }
   double dv(h*h*h), worst(0.0);
   for (int i = 0; i < nCheck; ++i) for (int j = 0; j <= i; ++j) {
       double expected(i == j ? 1.0 : 0.0);
       worst = std::max(worst, std::abs(S[i*nCheck+j]*dv - expected));
   }
   check(orbs->nAlpha() == nAlpha && orbs->nBeta() == nAlpha - (multiplicity - 1),
      "the numbers of alpha and beta electrons are those of the valence electrons");
   delete orbs;
   return worst;
}

int main(int argc, char** argv)
{
   QCoreApplication app(argc, argv);
   Calc::XtbProcess xtb(argc > 1 ? argv[1] : "");

   Calc::Input water;
   water.atomicNumbers << 8 << 1 << 1;
   water.coordinates << 0.9264 << -0.0675 << 0.0010
                     << 1.8854 << -0.0817 << -0.0104
                     << 0.6239 << -0.7774 << -0.5686;
   water.orbitals = true;
   Calc::Result r(xtb.energy(water));
   check(r.success && r.hasOrbitals, "water: orbitals are returned");
   if (r.hasOrbitals) {
      check(r.orbitals.nBasis() == 6 && r.orbitals.shells.size() == 4, "water: 6 basis functions");
      check(r.orbitals.coefficients.size() == 6, "water: 6 orbitals");
      check(r.orbitals.energies[0] < r.orbitals.energies[1] && r.orbitals.occupations[3] > 1.9
         && r.orbitals.occupations[4] < 0.1, "water: ordered by energy, four valence orbitals occupied");
      double dev(orthonormality(r, 6, 4));
      std::printf("   largest deviation from orthonormality: %.2e\n", dev);
      check(dev < 5e-3, "water: all orbitals are orthonormal (s and p functions)");
   }

   // sulfur has d functions
   Calc::Input h2s;
   h2s.atomicNumbers << 16 << 1 << 1;
   h2s.coordinates << 0.0 << 0.0 << 0.1 << 0.1 << 0.96 << -0.85 << -0.1 << -0.96 << -0.85;
   h2s.orbitals = true;
   Calc::Result s(xtb.energy(h2s));
   check(s.success && s.hasOrbitals && s.orbitals.nBasis() == 12, "H2S: 12 basis functions with d");
   if (s.hasOrbitals) {
      double dev(orthonormality(s, 8, 4));
      std::printf("   largest deviation from orthonormality: %.2e\n", dev);
      check(dev < 5e-3, "H2S: the orbitals are orthonormal (Cartesian d functions)");
   }

   // not requested, not returned
   water.orbitals = false;
   check(!xtb.energy(water).hasOrbitals, "no orbitals unless requested");

   std::printf("%d failure(s)\n", s_failures);
   return s_failures == 0 ? 0 : 1;
}
