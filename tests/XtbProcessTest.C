// Copyright (C) 2026 Roberto Peverati.  Part of IQmol, free software under the GNU General Public
// License, version 3 or (at your option) any later version.  Developed with the assistance of
// Claude (Anthropic).
// Standalone test for the xtb interface: optimizes water and checks the result.
// Usage: XtbProcessTest [path/to/xtb]     (or set XTB_PATH)

#include "Util/XtbProcess.h"
#include <QCoreApplication>
#include <QDir>
#include <QFile>
#include <QTemporaryDir>
#include <QRegularExpression>
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
#ifdef Q_OS_WIN
   QString fake(dir.filePath("xtb.bat"));      // Windows cannot start a shell script
   char const* sleeper("@echo off\r\nping -n 30 127.0.0.1 > nul\r\n");
#else
   QString fake(dir.filePath("xtb"));
   char const* sleeper("#!/bin/sh\nsleep 30\n");
#endif
   QFile script(fake);
   if (script.open(QIODevice::WriteOnly)) {
      script.write(sleeper);
      script.close();
      script.setPermissions(QFile::permissions(fake) | QFile::ExeOwner);
   }
   XtbProcess slow(fake);
   slow.setInteractive(true);
   QTimer::singleShot(300, [&slow]() { slow.cancel(); });
   Result cancelled(slow.optimize(water));
   check(cancelled.cancelled && !cancelled.success, "cancel() stops a running calculation");
   check(cancelled.seconds < 5.0, "cancellation returns promptly");

   // Vibrational frequencies of optimized water (values from the xtb program output)
   Input minimum(water);
   minimum.coordinates = opt.coordinates;
   Result freq(xtb.frequencies(minimum));
   check(freq.success && freq.hasVibrations && freq.vibrations.frequencies.size() == 3, "frequencies: three modes for water");
   if (freq.success && freq.vibrations.frequencies.size() == 3) {
      VibrationalData const& v(freq.vibrations);
      // (the exact values depend on how tightly the geometry was optimized, by a few cm-1)
      check(qAbs(v.frequencies[0] - 1539.0) < 5.0 && qAbs(v.frequencies[1] - 3645.0) < 5.0 
         && qAbs(v.frequencies[2] - 3653.0) < 5.0, "frequencies: water bend and stretches");
      check(v.intensities.size() == 3 && qAbs(v.intensities[0] - 133.0) < 3.0, "frequencies: IR intensity of the bend");
      check(qAbs(v.zeroPointEnergy - 0.02013) < 5e-5, "frequencies: zero point energy");
      double norm(0.0);
      for (double x : v.modes[0]) norm += x * x;
      check(qAbs(norm - 1.0) < 1e-9, "frequencies: normal modes are normalized");
      // mass weighted modes are orthogonal, which the Cartesian displacements are not
      double orth(0.0);
      static double const mass[3] = { 15.99940492, 1.00794075, 1.00794075 };
      for (int a = 0; a < 3; ++a) for (int b = a + 1; b < 3; ++b) {
          double sum(0.0);
          for (int i = 0; i < 9; ++i) sum += mass[i / 3] * v.modes[a][i] * v.modes[b][i];
          orth = qMax(orth, qAbs(sum));
      }
      check(orth < 2e-2, "frequencies: the modes are Cartesian displacements (orthogonal when mass weighted)");
   }
   check(!freq.success || freq.energy != 0.0, "frequencies: the energy is reported");

   // Implicit solvent (ALPB)
   Input inWater(water);
   inWater.solvent = "water";
   Result solvated(xtb.energy(inWater));
   check(solvated.success && solvated.hasGsolv && solvated.gsolv < 0.0, "ALPB(water) reports a negative Gsolv");
   check(solvated.success && solvated.energy < sp.energy, "water is stabilised by ALPB(water)");
   check(!sp.hasGsolv, "gas phase reports no Gsolv");

   Input inMeoh(water);
   inMeoh.solvent = "methanol";
   Result inMethanol(xtb.energy(inMeoh));
   check(inMethanol.success && inMethanol.energy != solvated.energy, "different solvents give different energies");

   Input inBad(water);
   inBad.solvent = "notasolvent";
   Result badSolvent(xtb.energy(inBad));
   check(!badSolvent.success && !badSolvent.message.isEmpty(), "unknown solvent is reported as failure");

   Result optWater(xtb.optimize(inWater));
   check(optWater.success && optWater.converged && optWater.hasGsolv, "optimization in ALPB(water) converges");

   // A temporary directory that is kept, for the xtb window to save its files
   {
      XtbProcess p(argc > 1 ? argv[1] : "");
      p.setKeepTemporaryDirectory(true);
      Result o(p.optimize(water));
      check(o.success && o.temporaryDirectory && QFile::exists(o.directory + "/xtbopt.xyz"), 
         "the temporary directory is kept and reported");
      Input minimum(water);
      minimum.coordinates = o.coordinates;
      Result h(p.frequencies(minimum));
      check(h.success && h.directory == o.directory && QFile::exists(o.directory + "/xtb_opt.out")
         && QFile::exists(o.directory + "/xtb_hess.out") && QFile::exists(o.directory + "/xtbopt.xyz")
         && QFile::exists(o.directory + "/g98.out"), "the optimization and the frequencies share the directory and keep their files");
      QDir(o.directory).removeRecursively();
      XtbProcess q(argc > 1 ? argv[1] : "");
      check(q.energy(water).directory.isEmpty(), "by default no directory is kept");
   }

   // Working directory and extra options
   {
      QTemporaryDir parent;
      QString work(parent.filePath("run/xtb files"));    // does not exist yet
      XtbProcess custom(argc > 1 ? argv[1] : "");
      custom.setWorkingDirectory(work);
      custom.setExtraArguments(QStringList() << "--etemp" << "5000" << "--chrg" << "1" << "--uhf" << "1");
      Result r(custom.optimize(water));
      check(r.success && r.converged, "working directory: optimization runs");
      check(QFile::exists(work + "/input.xyz") && QFile::exists(work + "/xtbopt.xyz"),
         "working directory: the input and output files are left there");
      QFile out(work + "/xtb_opt.out");
      QString text;
      if (out.open(QIODevice::ReadOnly)) text = QString::fromLocal8Bit(out.readAll());
      check(text.contains("TOTAL ENERGY"), "working directory: the output is saved as xtb_opt.out");
      check(text.contains("5000"), "extra options are passed to xtb (--etemp 5000)");
      check(r.success && r.energy > opt.energy + 0.1, "the extra options come last and override (water cation)");

      // a second run in the same directory must not pick up the first run's files
      QFile::copy(work + "/xtbopt.xyz", work + "/keep.xyz");
      Result again(custom.energy(water));
      check(again.success && !QFile::exists(work + "/xtbopt.xyz"), "stale results are removed before a run");
   }

   std::printf("%d failure(s)\n", s_failures);
   return s_failures == 0 ? 0 : 1;
}
