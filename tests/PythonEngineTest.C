// Copyright (C) 2026 Roberto Peverati.  Part of IQmol, free software under the GNU General Public
// License, version 3 or (at your option) any later version.  Developed with the assistance of
// Claude (Anthropic).
// Tests Calc::PythonEngine against the demo engine (a Morse toy potential) of the
// iqmolx-engines helper.
//
// Usage: PythonEngineTest <python> <python/iqmolx-engines/src>

#include "Util/PythonEngine.h"
#include <QCoreApplication>
#include <QElapsedTimer>
#include <QJsonArray>
#include <QTimer>
#include <QtMath>
#include <cstdio>

using namespace IQmol::Calc;

static int s_failures(0);

static void check(bool ok, QString const& what)
{
   std::printf("%s: %s\n", ok ? "PASS" : "FAIL", qPrintable(what));
   if (!ok) ++s_failures;
}

static double const eVPerHartree(27.211386245988);

static Input water(double stretch = 1.0)
{
   Input input;
   input.atomicNumbers << 8 << 1 << 1;
   input.coordinates << 0.0 << 0.0 << 0.117 
                     << 0.0 << 0.757*stretch << -0.469 
                     << 0.0 << -0.757*stretch << -0.469;
   return input;
}

static Input hydrogen(double r)
{
   Input input;
   input.atomicNumbers << 1 << 1;
   input.coordinates << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << r;
   return input;
}


int main(int argc, char** argv)
{
   QCoreApplication app(argc, argv);
   if (argc < 3) { std::printf("usage: PythonEngineTest <python> <helper src dir>\n"); return 2; }

   // A missing interpreter is reported, not fatal
   {
      PythonEngine engine("/nonexistent/python", "demo");
      Result r(engine.energy(hydrogen(0.62)));
      check(!r.success && r.message.contains("Unable to start"), "missing interpreter: " + r.message);
   }

   PythonEngine engine(argv[1], "demo");
   engine.setPythonPath(argv[2]);

   QString error;
   QJsonArray engines(engine.listEngines(&error));
   check(!engines.isEmpty(), "engines listed " + error);
   check(engine.engineInfo().value("name").toString() == "demo", "descriptor of the selected engine");
   check(engine.engineInfo().value("properties").toArray().size() == 5, "demo offers five properties");
   check(engine.name().startsWith("Demo"), "name from the title: " + engine.name());

   // H2 at the minimum of the Morse potential: E = -D = -4 eV, no force
   {
      Result r(engine.energy(hydrogen(0.62)));
      check(r.success, "H2 energy runs " + r.message);
      check(qAbs(r.energy * eVPerHartree + 4.0) < 1e-9, "H2 energy is -D");
      check(r.gradientNorm >= 0.0 && r.gradientNorm < 1e-9, "no gradient at the minimum");
      check(r.dipoleDebye >= 0.0, "dipole reported");
   }

   // Unsupported element is an error from the helper, which stays usable
   {
      Input uranium(hydrogen(0.62));
      uranium.atomicNumbers[0] = 92;
      Result r(engine.energy(uranium));
      check(!r.success && r.message.contains("Unsupported element"), "unsupported element: " + r.message);
      check(engine.energy(hydrogen(0.62)).success, "helper still works after an error");
   }

   // Optimization of a stretched water converges to a stationary point
   {
      Input input(water(1.4));
      Result r(engine.optimize(input));
      check(r.success && r.converged, "optimization converges " + r.message);
      check(r.energy < r.initialEnergy, "optimization lowers the energy");
      check(r.gradientNorm < 2e-3, QString("small final gradient %1").arg(r.gradientNorm));
      check(r.coordinates.size() == 9, "coordinates returned");
   }

   // Frequencies of H2 against the analytic Morse value  w = a sqrt(2 D / mu) / (2 pi c)
   {
      double const amu(1.66053906660e-27), eV(1.602176634e-19), c(2.99792458e10);
      double const mu(0.5 * 1.00782503223 * amu);
      double const expected(1.8e10 * qSqrt(2.0 * 4.0 * eV / mu) / (2.0 * M_PI * c));

      Result analytic(engine.frequencies(hydrogen(0.62)));
      check(analytic.success && analytic.hasVibrations, "frequencies (analytic Hessian) " + analytic.message);
      double top(analytic.vibrations.frequencies.isEmpty() ? 0.0 : analytic.vibrations.frequencies.last());
      check(qAbs(top / expected - 1.0) < 0.01, 
         QString("H2 frequency %1 cm-1, expected %2").arg(top).arg(expected));

      engine.setAnalyticHessian(false);
      Result numerical(engine.frequencies(hydrogen(0.62)));
      engine.setAnalyticHessian(true);
      check(numerical.success && numerical.hasVibrations, "frequencies (numerical Hessian) " + numerical.message);
      double top2(numerical.vibrations.frequencies.isEmpty() ? 0.0 : numerical.vibrations.frequencies.last());
      check(qAbs(top2 / top - 1.0) < 1e-3, 
         QString("numerical %1 and analytic %2 agree").arg(top2).arg(top));
   }

   // Cancelling a calculation that is waiting in the helper
   {
      QJsonObject slow;
      slow["delay_ms"] = 30000;
      engine.setComputeOptions(slow);
      engine.setInteractive(true);
      QTimer::singleShot(400, [&]() { engine.cancel(); });
      QElapsedTimer clock;
      clock.start();
      Result r(engine.energy(hydrogen(0.62)));
      check(r.cancelled && !r.success, "cancelled calculation is reported");
      check(clock.elapsed() < 10000, "cancel does not wait for the computation");

      engine.setComputeOptions(QJsonObject());
      Result again(engine.energy(hydrogen(0.62)));
      check(again.success, "helper usable after a cancel " + again.message);
      engine.setInteractive(false);
   }

   // Other load options load the model again
   {
      QJsonObject options;
      options["well_depth"] = 2.0;
      engine.setLoadOptions(options);
      Result r(engine.energy(hydrogen(0.62)));
      check(r.success && qAbs(r.energy * eVPerHartree + 2.0) < 1e-9, "load options take effect");
   }

   // Unknown engine
   {
      PythonEngine other(argv[1], "nonsense");
      other.setPythonPath(argv[2]);
      Result r(other.energy(hydrogen(0.62)));
      check(!r.success && r.message.contains("nonsense"), "unknown engine: " + r.message);
   }

   std::printf("\n%s\n", s_failures ? "FAILED" : "all tests passed");
   return s_failures ? 1 : 0;
}
