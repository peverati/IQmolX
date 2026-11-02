// Standalone test for the xtb interface: optimizes water and checks the result.
// Usage: XtbProcessTest [path/to/xtb]     (or set XTB_PATH)

#include "Util/XtbProcess.h"
#include <QCoreApplication>
#include <QFile>
#include <QTemporaryDir>
#include <QTimer>
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

   // Details parsed from the xtb output
   check(opt.success && opt.cycles > 0, "number of optimization cycles is parsed");
   check(opt.success && opt.gradientNorm >= 0.0 && opt.gradientNorm < 0.01, "gradient norm is parsed");
   check(opt.success && qAbs(opt.gapEV - 14.4) < 0.5, "HOMO-LUMO gap is parsed");
   check(opt.success && qAbs(opt.dipoleDebye - 2.2) < 0.3, "dipole moment is parsed");
   check(opt.success && opt.initialEnergy > opt.energy, "initial energy is parsed and above the final one");
   check(opt.success && opt.seconds > 0.0, "wall time is recorded");

   // Interactive mode: event loop running, output delivered as it arrives
   QString streamed;
   XtbProcess live(argc > 1 ? argv[1] : "");
   live.setInteractive(true);
   live.setProgressCallback([&streamed](QString const& chunk) { streamed += chunk; });
   Result liveResult(live.optimize(water));
   check(liveResult.success && liveResult.converged, "interactive optimization converges");
   check(streamed.contains("GEOMETRY OPTIMIZATION CONVERGED"), "progress callback receives the output");
   check(qAbs(liveResult.energy - opt.energy) < 1e-6, "interactive and blocking modes agree");

   // Cancellation, using a stand-in for xtb that just sleeps
   QTemporaryDir dir;
   QString fake(dir.filePath("xtb"));
   QFile script(fake);
   if (script.open(QIODevice::WriteOnly)) {
      script.write("#!/bin/sh\nsleep 30\n");
      script.close();
      script.setPermissions(QFile::permissions(fake) | QFile::ExeOwner);
   }
   XtbProcess slow(fake);
   slow.setInteractive(true);
   QTimer::singleShot(300, [&slow]() { slow.cancel(); });
   Result cancelled(slow.optimize(water));
   check(cancelled.cancelled && !cancelled.success, "cancel() stops a running calculation");
   check(cancelled.seconds < 5.0, "cancellation returns promptly");

   std::printf("%d failure(s)\n", s_failures);
   return s_failures == 0 ? 0 : 1;
}
