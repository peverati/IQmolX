// Compares the two xtb backends, the external program (XtbProcess) and libxtb with IQmol's
// own optimizer (XtbLibrary), on a set of molecules in tests/data.
//
// Usage: XtbBackendTest <path/to/xtb> <path/to/tests/data>

#include "Util/XtbLibrary.h"
#include "Util/XtbProcess.h"
#include <QFile>
#include <QStringList>
#include <QTextStream>
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

   std::printf("%d failure(s)\n", s_failures);
   return s_failures == 0 ? 0 : 1;
}
