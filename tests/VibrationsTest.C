// Tests of the harmonic vibrational analysis on cases with analytic answers.

#include "Util/Vibrations.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>

using namespace IQmol::Calc;
typedef std::vector<double> Vec;

static int s_failures(0);

static void check(bool ok, char const* what)
{
   std::printf("%s: %s\n", ok ? "PASS" : "FAIL", what);
   if (!ok) ++s_failures;
}

static double const s_factor(5140.487446);

// Hessian of E = 1/2 sum_pairs k (r - r0)^2 at equilibrium: each pair adds k u u^T
static Vec pairSpringHessian(Vec const& x, double k, bool allPairs = true, int skipI = -1, int skipJ = -1)
{
   int n(x.size() / 3), dim(3 * n);
   Vec h(dim * dim, 0.0);
   for (int i = 0; i < n; ++i) {
       for (int j = i + 1; j < n; ++j) {
           if (!allPairs) continue;
           if (i == skipI && j == skipJ) continue;
           double u[3], r(0.0);
           for (int c = 0; c < 3; ++c) { u[c] = x[3*j+c] - x[3*i+c]; r += u[c] * u[c]; }
           r = std::sqrt(r);
           for (int c = 0; c < 3; ++c) u[c] /= r;
           for (int a = 0; a < 3; ++a) {
               for (int b = 0; b < 3; ++b) {
                   double v(k * u[a] * u[b]);
                   h[(3*i+a)*dim + 3*i+b] += v;  h[(3*j+a)*dim + 3*j+b] += v;
                   h[(3*i+a)*dim + 3*j+b] -= v;  h[(3*j+a)*dim + 3*i+b] -= v;
               }
           }
       }
   }
   return h;
}

// Rotates the coordinates and the Hessian (block-wise) by the same rotation
static void rotate(Vec& x, Vec& h, double angleZ, double angleX)
{
   double cz(std::cos(angleZ)), sz(std::sin(angleZ)), cx(std::cos(angleX)), sx(std::sin(angleX));
   double Rz[3][3] = { {cz, -sz, 0}, {sz, cz, 0}, {0, 0, 1} };
   double Rx[3][3] = { {1, 0, 0}, {0, cx, -sx}, {0, sx, cx} };
   double R[3][3];
   for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
      R[a][b] = 0; for (int k = 0; k < 3; ++k) R[a][b] += Rx[a][k] * Rz[k][b];
   }
   int n(x.size() / 3), dim(3 * n);
   for (int i = 0; i < n; ++i) {
       double v[3] = { x[3*i], x[3*i+1], x[3*i+2] };
       for (int a = 0; a < 3; ++a) { x[3*i+a] = 0; for (int b = 0; b < 3; ++b) x[3*i+a] += R[a][b] * v[b]; }
   }
   Vec out(h.size(), 0.0);
   for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j)
      for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
         double sum(0.0);
         for (int c = 0; c < 3; ++c) for (int d = 0; d < 3; ++d)
            sum += R[a][c] * h[(3*i+c)*dim + 3*j+d] * R[b][d];
         out[(3*i+a)*dim + 3*j+b] = sum;
      }
   h = out;
}


int main()
{
   // ---- masses
   check(std::fabs(atomicMass(1) - 1.00794075) < 1e-8 && std::fabs(atomicMass(8) - 15.99940492) < 1e-8,
      "the atomic masses are those of xtb");
   check(atomicMass(0) == 0.0 && atomicMass(200) == 0.0, "unknown elements have no mass");

   // ---- eigenvalue solver
   {
      int n(25);
      Vec a(n * n), original;
      srand(7);
      for (int i = 0; i < n; ++i) for (int j = i; j < n; ++j) a[i*n+j] = a[j*n+i] = (rand() / double(RAND_MAX)) - 0.5;
      original = a;
      Vec values, vectors;
      symmetricEigenvalues(a, n, values, vectors);
      double residual(0.0), orthogonality(0.0);
      for (int k = 0; k < n; ++k) {
          for (int i = 0; i < n; ++i) {
              double sum(0.0);
              for (int j = 0; j < n; ++j) sum += original[i*n+j] * vectors[j*n+k];
              residual = std::fmax(residual, std::fabs(sum - values[k] * vectors[i*n+k]));
          }
          for (int l = 0; l < n; ++l) {
              double sum(0.0);
              for (int i = 0; i < n; ++i) sum += vectors[i*n+k] * vectors[i*n+l];
              orthogonality = std::fmax(orthogonality, std::fabs(sum - (k == l ? 1.0 : 0.0)));
          }
      }
      bool sorted(true);
      for (int k = 1; k < n; ++k) if (values[k] < values[k-1]) sorted = false;
      std::printf("  eigen solver: residual %.1e, orthogonality %.1e\n", residual, orthogonality);
      check(residual < 1e-12 && orthogonality < 1e-12 && sorted, "eigenvalues and vectors of a symmetric matrix");
   }

   // ---- diatomic: nu = factor * sqrt(k / mu)
   {
      double k(0.37), m1(atomicMass(1)), m2(atomicMass(9));
      Vec x = { 0.0, 0.0, 0.0,   0.0, 0.0, 1.7 };
      Vec h(36, 0.0);
      for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) h[(3*a+2)*6 + 3*b+2] = (a == b) ? k : -k;
      std::vector<int> Z = { 1, 9 };
      VibrationalAnalysis v(analyseVibrations(Z, x, h));
      double mu(m1 * m2 / (m1 + m2));
      std::printf("  diatomic: %.6f cm-1, expected %.6f\n", v.frequencies.empty() ? 0.0 : v.frequencies[0], s_factor * std::sqrt(k / mu));
      check(v.ok && v.frequencies.size() == 1 && std::fabs(v.frequencies[0] - s_factor * std::sqrt(k / mu)) < 1e-6,
         "diatomic frequency is sqrt(k/mu)");
      check(v.translationsRotations == 5, "a diatomic has 5 translations and rotations");
      // The reduced mass is defined as in Gaussian and xtb, 1/sum(l^2) with l the Cartesian
      // displacements normalized in mass weighted coordinates (H2 has 1.008, not 0.504)
      double gaussianMass((1.0 / m1 + 1.0 / m2) / (1.0 / (m1 * m1) + 1.0 / (m2 * m2)));
      check(v.ok && std::fabs(v.reducedMasses[0] - gaussianMass) < 1e-9, "reduced mass in the Gaussian convention");
      // the displacements keep the centre of mass fixed and are normalized
      double cm(m1 * v.modes[0][2] + m2 * v.modes[0][5]), norm(0.0);
      for (double d : v.modes[0]) norm += d * d;
      check(std::fabs(cm) < 1e-10 && std::fabs(norm - 1.0) < 1e-12, "the mode is normalized and does not move the centre of mass");

      // IR intensity: dmu/dQ = a / sqrt(mu) for dmu_z/dz = (+a, -a)
      double a(0.5);
      Vec dip(18, 0.0);
      dip[3*2 + 2] = a;  dip[3*5 + 2] = -a;
      VibrationalAnalysis w(analyseVibrations(Z, x, h, &dip));
      double expected(974.8801 * a * a / mu);
      std::printf("  intensity: %.4f km/mol, expected %.4f\n", w.intensities.empty() ? 0.0 : w.intensities[0], expected);
      check(w.ok && w.intensities.size() == 1 && std::fabs(w.intensities[0] - expected) < 1e-6 * expected,
         "IR intensity of a diatomic");
      check(std::fabs(w.zeroPointEnergy - 0.5 * v.frequencies[0] * 4.556335252912e-6) < 1e-12, "zero point energy");
   }

   // ---- bent triatomic with springs between all pairs
   Vec bentX, bentH;
   std::vector<int> bentZ = { 8, 1, 1 };
   {
      bentX = { 0.0, 0.0, 0.0,   0.0, 1.43, 1.10,   0.0, -1.43, 1.10 };
      bentH = pairSpringHessian(bentX, 0.4);
      VibrationalAnalysis v(analyseVibrations(bentZ, bentX, bentH));
      check(v.ok && v.frequencies.size() == 3 && v.translationsRotations == 6, "a bent triatomic has 3 vibrations");
      bool real(true);
      for (double f : v.frequencies) if (f <= 0.0) real = false;
      check(real, "a harmonic minimum has no imaginary frequencies");

      // Rotating the molecule does not change the frequencies
      Vec x(bentX), h(bentH);
      rotate(x, h, 0.7, 1.9);
      VibrationalAnalysis r(analyseVibrations(bentZ, x, h));
      double diff(0.0);
      for (int i = 0; i < 3; ++i) diff = std::fmax(diff, std::fabs(r.frequencies[i] - v.frequencies[i]));
      std::printf("  rotation changes the frequencies by %.1e cm-1\n", diff);
      check(r.ok && diff < 1e-8, "frequencies do not depend on the orientation");

      // Translating the molecule does not change them either
      Vec shifted(bentX);
      for (size_t i = 0; i < shifted.size(); ++i) shifted[i] += (i % 3 == 0 ? 3.0 : (i % 3 == 1 ? -2.0 : 5.0));
      VibrationalAnalysis s(analyseVibrations(bentZ, shifted, bentH));
      double diff2(0.0);
      for (int i = 0; i < 3; ++i) diff2 = std::fmax(diff2, std::fabs(s.frequencies[i] - v.frequencies[i]));
      check(s.ok && diff2 < 1e-8, "frequencies do not depend on the position");

      // Mass weighted modes are orthogonal: sum m d_a d_b = 0
      double orth(0.0);
      for (int a = 0; a < 3; ++a) for (int b = a + 1; b < 3; ++b) {
          double sum(0.0);
          for (int i = 0; i < 9; ++i) sum += atomicMass(bentZ[i/3]) * v.modes[a][i] * v.modes[b][i];
          orth = std::fmax(orth, std::fabs(sum));
      }
      check(orth < 1e-10, "the normal modes are orthogonal in mass weighted coordinates");

      // A negative curvature appears as an imaginary (negative) frequency
      Vec soft(bentH);
      Vec hb(pairSpringHessian(bentX, -0.4, true, 1, 2));
      for (size_t i = 0; i < soft.size(); ++i) soft[i] = bentH[i] + 2.0 * (hb[i]);
      VibrationalAnalysis im(analyseVibrations(bentZ, bentX, soft));
      bool anyNegative(false);
      for (double f : im.frequencies) if (f < 0.0) anyNegative = true;
      check(im.ok && anyNegative && im.frequencies[0] < 0.0, "negative curvature gives imaginary frequencies first");
      check(im.ok && im.zeroPointEnergy < v.zeroPointEnergy + 1.0, "imaginary modes do not enter the zero point energy");
   }

   // ---- linear triatomic: 3N-5 modes, and the zero bends are not mistaken for rotations
   {
      Vec x = { -2.2, 0.0, 0.0,   0.0, 0.0, 0.0,   2.2, 0.0, 0.0 };
      Vec h(pairSpringHessian(x, 0.5, true, 0, 2));      // no spring between the ends
      std::vector<int> Z = { 8, 6, 8 };
      VibrationalAnalysis v(analyseVibrations(Z, x, h));
      int zero(0);
      for (double f : v.frequencies) if (std::fabs(f) < 1e-3) ++zero;
      check(v.ok && v.frequencies.size() == 4 && v.translationsRotations == 5, "a linear triatomic has 3N-5 = 4 modes");
      check(zero == 2, "the two bending modes of the springs-only model have zero frequency");
   }

   // ---- single atom, bad input
   {
      Vec x = { 0, 0, 0 };
      Vec h(9, 0.0);
      std::vector<int> Z = { 6 };
      VibrationalAnalysis v(analyseVibrations(Z, x, h));
      check(v.ok && v.frequencies.empty(), "a single atom has no vibrations");

      std::vector<int> Z2 = { 1, 1 };
      VibrationalAnalysis bad(analyseVibrations(Z2, x, h));
      check(!bad.ok && !bad.message.empty(), "inconsistent sizes are reported");

      std::vector<int> Z3 = { 0 };
      VibrationalAnalysis unknown(analyseVibrations(Z3, x, h));
      check(!unknown.ok && !unknown.message.empty(), "an unknown element is reported");
   }

   std::printf("%d failure(s)\n", s_failures);
   return s_failures == 0 ? 0 : 1;
}
