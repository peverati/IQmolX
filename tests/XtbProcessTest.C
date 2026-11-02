// Standalone test for the xtb interface: optimizes water and checks the result.
// Usage: XtbProcessTest [path/to/xtb]     (or set XTB_PATH)

#include "Util/XtbProcess.h"
#include <QCoreApplication>
#include <QtMath>
#include <cstdio>

using namespace IQmol::Calc;

static int s_failures(0);

static void check(bool ok, char const* what)
{
   std::printf("%s: %s\n", ok ? "PASS" : "FAIL", what);
   if (!ok) ++s_failures;
}

static double distance(QVector<double> const& c, int i, int j)
{
   double dx(c[3*i]-c[3*j]), dy(c[3*i+1]-c[3*j+1]), dz(c[3*i+2]-c[3*j+2]);
   return qSqrt(dx*dx + dy*dy + dz*dz);
}


int main(int argc, char** argv)
{
   QCoreApplication app(argc, argv);

   XtbProcess xtb(argc > 1 ? argv[1] : "");

   Input water;
   water.atomicNumbers << 8 << 1 << 1;
   water.coordinates << 0.0 <<  0.0 << 0.0
                     << 0.0 <<  0.9 << 0.6
                     << 0.0 << -0.9 << 0.6;

   Result sp(xtb.energy(water));
   std::printf("single point: %.8f Eh  %s\n", sp.energy, qPrintable(sp.message));
   check(sp.success, "single point runs");
   check(sp.success && qAbs(sp.energy + 5.0705) < 0.05, "single point energy is sensible");

   Result opt(xtb.optimize(water));
   std::printf("optimized:    %.8f Eh  %s\n", opt.energy, qPrintable(opt.message));
   check(opt.success && opt.converged, "optimization converges");
   if (opt.success) {
      double r(distance(opt.coordinates, 0, 1));
      std::printf("O-H = %.4f A\n", r);
      check(r > 0.94 && r < 1.00, "O-H distance is near 0.97 A");
      check(opt.energy <= sp.energy + 1e-9, "optimization lowers the energy");
   }

   // Charged and open-shell: H2O+ doublet
   Input cation(water);
   cation.charge = 1;
   cation.multiplicity = 2;
   Result rad(xtb.energy(cation));
   check(rad.success && rad.energy > sp.energy, "cation doublet runs, higher in energy");

   // Bad input is reported, not crashed
   Input bad(water);
   bad.multiplicity = 2;   // 10 electrons cannot be a doublet
   Result fail(xtb.energy(bad));
   check(!fail.success && !fail.message.isEmpty(), "inconsistent spin is reported as failure");

   std::printf("%d failure(s)\n", s_failures);
   return s_failures == 0 ? 0 : 1;
}
