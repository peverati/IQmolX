// Copyright (C) 2026 Roberto Peverati.  Part of IQmol, free software under the GNU General Public
// License, version 3 or (at your option) any later version.
// Compares the two xtb backends, the external program (XtbProcess) and libxtb with IQmol's
// own optimizer (XtbLibrary), on a set of molecules in tests/data.
//
// Usage: XtbBackendTest <path/to/xtb> <path/to/tests/data>

#include "Util/XtbLibrary.h"
#include "Util/XtbProcess.h"
#include <QFile>
#include <QStringList>
#include <QTextStream>
#include <QElapsedTimer>
#include <QtMath>
#include <cstdio>

using namespace IQmol::Calc;

static int s_failures(0);

static void check(bool ok, QString const& what)
{
   std::printf("%s: %s\n", ok ? "PASS" : "FAIL", qPrintable(what));
   if (!ok) ++s_failures;
}

static int atomicNumber(QString const& symbol)
{
   static QStringList symbols = QString("X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca").split(' ');
   return symbols.indexOf(symbol);
}

static bool readXyz(QString const& path, Input& input)
{
   QFile file(path);
   if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return false;
   QTextStream in(&file);
   int n(in.readLine().trimmed().toInt());
   in.readLine();
   for (int i = 0; i < n; ++i) {
       QStringList t(in.readLine().split(' ', Qt::SkipEmptyParts));
       if (t.size() < 4) return false;
       input.atomicNumbers.append(atomicNumber(t[0]));
       input.coordinates << t[1].toDouble() << t[2].toDouble() << t[3].toDouble();
   }
   return true;
}

// Root-mean-square difference of all interatomic distances: independent of orientation
static double distanceRmsd(QVector<double> const& a, QVector<double> const& b)
{
   int n(a.size() / 3);
   double sum(0.0);
   int count(0);
   for (int i = 0; i < n; ++i) {
       for (int j = i + 1; j < n; ++j) {
           double da(0.0), db(0.0);
           for (int k = 0; k < 3; ++k) {
               da += qPow(a[3*i+k] - a[3*j+k], 2);
               db += qPow(b[3*i+k] - b[3*j+k], 2);
           }
           sum += qPow(qSqrt(da) - qSqrt(db), 2);
           ++count;
       }
   }
   return qSqrt(sum / count);
}


struct Case { char const* file; int charge; int multiplicity; char const* solvent; };


int main(int argc, char** argv)
{
   if (argc < 3) { std::printf("usage: XtbBackendTest <xtb> <data dir>\n"); return 2; }
   QString xtbPath(argv[1]), dataDir(argv[2]);

   Case cases[] = {
      { "water",              0, 1, ""      },
      { "ethanol",            0, 1, ""      },
      { "ammonium",           1, 1, ""      },
      { "methyl",             0, 2, ""      },
      { "butane",             0, 1, ""      },
      { "aspirin",            0, 1, "water" },
      { "caffeine",           0, 1, ""      },
      { "glycine_zwitterion", 0, 1, "water" },
   };

   std::printf("%-20s %5s  %14s %14s  %9s  %7s  %s\n", "molecule", "atoms", "E program",
      "E library", "dE (Eh)", "rmsd(A)", "cycles prog/lib");

   for (Case const& c : cases) {
       Input input;
       if (!readXyz(dataDir + "/" + c.file + ".xyz", input)) {
          check(false, QString("read ") + c.file);
          continue;
       }
       input.charge = c.charge;
       input.multiplicity = c.multiplicity;
       input.solvent = c.solvent;
       QString label(QString(c.file) + (input.solvent.isEmpty() ? "" : " in " + input.solvent));

       XtbProcess program(xtbPath);
       XtbLibrary library;

       // Single points: the same energy
       Result spP(program.energy(input)), spL(library.energy(input));
       check(spP.success && spL.success && qAbs(spP.energy - spL.energy) < 1e-6,
          label + ": single point energies agree (" + QString::number(spL.energy - spP.energy, 'g', 2) + " Eh)");
       if (!input.solvent.isEmpty()) {
          check(spP.hasGsolv && spL.hasGsolv && qAbs(spP.gsolv - spL.gsolv) < 1e-6, label + ": Gsolv agrees");
       }
       if (spP.success && spL.success) {
          check(qAbs(spP.gapEV - spL.gapEV) < 1e-3, label + ": HOMO-LUMO gap agrees");
          if (spP.dipoleDebye >= 0.0) {
             check(qAbs(spP.dipoleDebye - spL.dipoleDebye) < 2e-3, label + ": dipole agrees");
          }
       }

       // Optimizations: the same minimum
       Result optP(program.optimize(input)), optL(library.optimize(input));
       if (!optP.success || !optL.success) {
          check(false, label + ": both optimizations run: " + optP.message + " / " + optL.message);
          continue;
       }
       double dE(optL.energy - optP.energy);
       double rmsd(distanceRmsd(optP.coordinates, optL.coordinates));
       std::printf("%-20s %5d  %14.8f %14.8f  %9.2e  %7.4f  %d / %d\n", qPrintable(label),
          input.atomicNumbers.size(), optP.energy, optL.energy, dE, rmsd, optP.cycles, optL.cycles);
       check(optL.converged, label + ": library optimization converges");
       check(qAbs(dE) < 5e-5, label + ": same minimum energy");
       check(rmsd < 0.05, label + ": same geometry");
       check(optL.gradientNorm >= 0.0 && optL.gradientNorm < 1.5e-3, label + ": final gradient is small");

       // Vibrational frequencies of the same geometry from both backends
       Input minimum(input);
       minimum.coordinates = optP.coordinates;
       QElapsedTimer clockP;  clockP.start();
       Result freqP(program.frequencies(minimum));
       double secondsP(clockP.elapsed() / 1000.0);
       QElapsedTimer clockL;  clockL.start();
       Result freqL(library.frequencies(minimum));
       double secondsL(clockL.elapsed() / 1000.0);
       if (!freqP.success || !freqL.success) {
          check(false, label + ": both frequency calculations run: " + freqP.message + " / " + freqL.message);
          continue;
       }

       VibrationalData const& vp(freqP.vibrations);
       VibrationalData const& vl(freqL.vibrations);
       int natoms(input.atomicNumbers.size());
       check(vp.frequencies.size() == vl.frequencies.size() && vl.frequencies.size() >= 3 * natoms - 6,
          label + ": same number of modes (" + QString::number(vl.frequencies.size()) + ")");

       // The two Hessians differ slightly (the library's agrees with finite differences of
       // the gradient, the program's --hess deviates from them by about 0.05% in an O-H
       // stretch), so the comparison allows 0.2% plus 4 cm-1 for soft modes.  Modes that
       // lie close to a neighbour mix, so their shapes are only compared if they are
       // separated by more than 30 cm-1.
       double maxDf(0.0), maxRelDf(0.0), maxDi(0.0), maxDmass(0.0), minOverlap(1.0);
       double sumP(0.0), sumL(0.0);
       bool frequenciesAgree(true), intensitiesAgree(true);
       int count(qMin(vp.frequencies.size(), vl.frequencies.size())), compared(0);
       for (int i = 0; i < count; ++i) {
           double df(qAbs(vp.frequencies[i] - vl.frequencies[i]));
           maxDf = qMax(maxDf, df);
           maxRelDf = qMax(maxRelDf, df / qMax(100.0, qAbs(vp.frequencies[i])));
           if (df > 0.002 * qAbs(vp.frequencies[i]) + 4.0) frequenciesAgree = false;

           bool isolated(true);
           if (i > 0 && qAbs(vp.frequencies[i] - vp.frequencies[i-1]) < 30.0) isolated = false;
           if (i + 1 < count && qAbs(vp.frequencies[i] - vp.frequencies[i+1]) < 30.0) isolated = false;

           if (vp.intensities.size() == count && vl.intensities.size() == count) {
              sumP += vp.intensities[i];
              sumL += vl.intensities[i];
              if (isolated) {
                 double dI(qAbs(vp.intensities[i] - vl.intensities[i]));
                 maxDi = qMax(maxDi, dI);
                 if (dI > 0.05 * vp.intensities[i] + 0.6) intensitiesAgree = false;
              }
           }
           if (isolated) {
              maxDmass = qMax(maxDmass, qAbs(vp.reducedMasses[i] - vl.reducedMasses[i]) / vp.reducedMasses[i]);
           }
           if (isolated && vp.modes[i].size() == vl.modes[i].size()) {
              double overlap(0.0);
              for (int k = 0; k < vp.modes[i].size(); ++k) overlap += vp.modes[i][k] * vl.modes[i][k];
              minOverlap = qMin(minOverlap, qAbs(overlap));
              ++compared;
           }
       }
       std::printf("   frequencies: %d modes, %.1f to %.1f cm-1, max df %.2f cm-1 (%.3f%%), max dI %.2f km/mol, total IR %.1f/%.1f km/mol, "
                   "%d isolated modes compared, overlap >= %.4f, ZPE %.6f/%.6f Eh, imaginary %d/%d, "
                   "time %.2f s program / %.2f s library\n", count, vl.frequencies.first(),
                   vl.frequencies.last(), maxDf, 100.0 * maxRelDf, maxDi, sumP, sumL, compared, minOverlap, vp.zeroPointEnergy,
                   vl.zeroPointEnergy, vp.imaginary(), vl.imaginary(), secondsP, secondsL);
       check(frequenciesAgree, label + ": frequencies agree (0.2% + 4 cm-1)");
       check(intensitiesAgree, label + ": IR intensities of isolated modes agree (5% + 0.6 km/mol)");
       check(qAbs(sumP - sumL) < 0.03 * sumP + 1.0, label + ": total IR intensity agrees within 3%");
       check(maxDmass < 0.02, label + ": reduced masses of isolated modes agree within 2%");
       check(minOverlap > 0.97, label + ": isolated normal modes agree");
       check(qAbs(vp.zeroPointEnergy - vl.zeroPointEnergy) < 2e-4 * vp.zeroPointEnergy + 1e-5,
          label + ": zero point energies agree");
       check(vp.imaginary() == vl.imaginary(), label + ": both find the same number of imaginary frequencies");
       check(freqL.gradientNorm >= 0.0 && freqL.gradientNorm < 2e-3, label + ": the geometry is stationary");
   }

   // GFN-FF is a force field: it does not use the electrons, so any spin is accepted, 
   // and the electronic temperature is ignored
   {
      Input methyl;
      readXyz(dataDir + "/methyl.xyz", methyl);
      XtbProcess program(xtbPath, XtbProcess::GFNFF);
      XtbLibrary library(XtbProcess::GFNFF);
      for (int multiplicity : { 1, 2 }) {
         methyl.multiplicity = multiplicity;
         Result p(program.energy(methyl)), l(library.energy(methyl));
         check(p.success && l.success && qAbs(p.energy - l.energy) < 1e-6,
            QString("GFN-FF energies agree for a radical with multiplicity %1").arg(multiplicity));
      }
      methyl.electronicTemperature = 1000.0;
      Result hot(library.energy(methyl));
      methyl.electronicTemperature = -1.0;
      Result cold(library.energy(methyl));
      check(hot.success && cold.success && qAbs(hot.energy - cold.energy) < 1e-10, "GFN-FF ignores the electronic temperature");

      Input ammonium;
      readXyz(dataDir + "/ammonium.xyz", ammonium);
      ammonium.charge = 1;
      Result charged(library.optimize(ammonium)), chargedP(program.optimize(ammonium));
      check(charged.success && chargedP.success && qAbs(charged.energy - chargedP.energy) < 5e-5, "GFN-FF optimizations of a cation agree");
      Input overlap;
      overlap.atomicNumbers << 8 << 1 << 1;
      overlap.coordinates << 0.0 << 0.0 << 0.0 << 0.5 << 0.0 << 0.0 << 0.5 << 0.0 << 0.0;
      Result bad(library.energy(overlap));
      check(!bad.success && !bad.message.isEmpty(), "GFN-FF reports atoms on top of each other instead of crashing");
   }

   // The electronic temperature of the xTB methods
   {
      Input water;
      readXyz(dataDir + "/water.xyz", water);
      XtbProcess program(xtbPath);
      XtbLibrary library;
      Result standard(library.energy(water));                // the program default, 300 K
      water.electronicTemperature = 300.0;
      Result explicit300(library.energy(water)), explicit300P(program.energy(water));
      check(standard.success && qAbs(standard.energy - explicit300.energy) < 1e-10 && qAbs(explicit300P.energy - explicit300.energy) < 1e-6,
         "300 K is the default, in the library and in the program");
      water.electronicTemperature = 0.0;
      Result zero(library.energy(water)), zeroP(program.energy(water));
      check(zero.success && zeroP.success && qAbs(zero.energy - standard.energy) < 1e-8, "a closed-shell molecule is the same at 0 K");

      // A (nearly) degenerate HOMO and LUMO needs smearing: this must fail cleanly, not crash
      Input twisted;
      readXyz(dataDir + "/twisted_ethylene.xyz", twisted);
      twisted.electronicTemperature = 300.0;
      Result warm(library.energy(twisted)), warmP(program.energy(twisted));
      check(warm.success && warmP.success && qAbs(warm.energy - warmP.energy) < 1e-6, "twisted ethylene works at 300 K with both backends");
      twisted.electronicTemperature = 100.0;
      Result cooler(library.energy(twisted));
      check(cooler.success && qAbs(cooler.energy - warm.energy) > 1e-4, "the electronic temperature changes the energy of a small gap molecule");
      twisted.electronicTemperature = 0.0;
      Result off(library.energy(twisted)), offP(program.energy(twisted));
      check(!off.success && !off.message.isEmpty(), "at 0 K the library reports the failure instead of crashing: " + off.message.left(60));
      check(!offP.success && !offP.message.isEmpty(), "at 0 K the program reports the failure");
   }

   // Error handling and the features of the shared interface
   {
      Input water;
      readXyz(dataDir + "/water.xyz", water);
      XtbLibrary library;

      Input badSpin(water);
      badSpin.multiplicity = 2;
      Result r(library.energy(badSpin));
      check(!r.success && r.message.contains("inconsistent"), "library rejects an impossible spin: " + r.message);

      Input badSolvent(water);
      badSolvent.solvent = "notasolvent";
      Result s(library.energy(badSolvent));
      check(!s.success && !s.message.isEmpty(), "library reports an unknown solvent");

      QString streamed;
      library.setProgressCallback([&streamed](QString const& t) { streamed += t; });
      Result p(library.optimize(water));
      check(p.success && streamed.contains("cycle") && streamed.contains("CONVERGED"), "progress is reported");

      XtbLibrary stopper;
      stopper.setInteractive(true);
      int cycles(0);
      stopper.setProgressCallback([&](QString const& t) { if (t.contains("cycle") && ++cycles == 2) stopper.cancel(); });
      Input aspirin;
      readXyz(dataDir + "/aspirin.xyz", aspirin);
      Result c(stopper.optimize(aspirin));
      check(c.cancelled && !c.success, "cancel() stops the optimization");
   }

   // Orbitals: libxtb does not report the basis, so the library backend gets them from the program
   {
      Input water;
      readXyz(dataDir + "/water.xyz", water);
      water.orbitals = true;
      XtbLibrary library;
      Result none(library.optimize(water));
      check(none.success && !none.hasOrbitals, "library without the program: energy, no orbitals");

      library.setExecutable(xtbPath);
      Result r(library.optimize(water));
      check(r.success && r.hasOrbitals && r.orbitals.nBasis() == 6, "library backend returns orbitals via the program");

      XtbProcess program(xtbPath);
      Result reference(program.energy([&]() { Input i(water); i.coordinates = r.coordinates; return i; }()));
      check(r.hasOrbitals && reference.hasOrbitals 
         && qAbs(r.orbitals.energies[0] - reference.orbitals.energies[0]) < 1e-9,
         "the orbitals belong to the optimized geometry");
   }

   std::printf("%d failure(s)\n", s_failures);
   return s_failures == 0 ? 0 : 1;
}
