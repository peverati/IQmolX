// Tests of the L-BFGS optimizer on analytic functions (no xtb needed).

#include "Util/Optimizer.h"
#include <cstdio>

using namespace IQmol::Calc;

static int s_failures(0);

static void check(bool ok, char const* what)
{
   std::printf("%s: %s\n", ok ? "PASS" : "FAIL", what);
   if (!ok) ++s_failures;
}

typedef std::vector<double> Vec;


// A stiff/soft anisotropic quadratic with its minimum at (1,2,3,...)
static bool quadratic(Vec const& x, double& e, Vec& g)
{
   g.assign(x.size(), 0.0);
   e = 0.0;
   for (size_t i = 0; i < x.size(); ++i) {
       double k(0.05 + 2.0 * i);                 // curvatures from soft to stiff
       double d(x[i] - (1.0 + i));
       e += 0.5 * k * d * d;
       g[i] = k * d;
   }
   return true;
}


// Rosenbrock banana, minimum at (1,1)
static bool rosenbrock(Vec const& x, double& e, Vec& g)
{
   double a(1.0 - x[0]), b(x[1] - x[0] * x[0]);
   e = a * a + 100.0 * b * b;
   g.assign(2, 0.0);
   g[0] = -2.0 * a - 400.0 * x[0] * b;
   g[1] = 200.0 * b;
   return true;
}


// Lennard-Jones cluster (reduced units): 7 atoms, global minimum is the pentagonal bipyramid, E = -16.505384
static bool lennardJones(Vec const& x, double& e, Vec& g)
{
   size_t n(x.size() / 3);
   e = 0.0;
   g.assign(x.size(), 0.0);
   for (size_t i = 0; i < n; ++i) {
       for (size_t j = i + 1; j < n; ++j) {
           double d[3], r2(0.0);
           for (int k = 0; k < 3; ++k) { d[k] = x[3*i+k] - x[3*j+k]; r2 += d[k] * d[k]; }
           double inv6(1.0 / (r2 * r2 * r2));
           e += 4.0 * (inv6 * inv6 - inv6);
           double f(4.0 * (-12.0 * inv6 * inv6 + 6.0 * inv6) / r2);   // dE/dr * 1/r
           for (int k = 0; k < 3; ++k) { g[3*i+k] += f * d[k]; g[3*j+k] -= f * d[k]; }
       }
   }
   return true;
}


int main()
{
   OptimizerSettings tight;
   tight.gradientTol = 1e-6;
   tight.energyTol = 1e-12;
   tight.maxCycles = 2000;

   {
      Vec x0(6, 0.0);
      OptimizerResult r(Optimizer(tight).minimize(quadratic, x0));
      double error(0.0);
      for (size_t i = 0; i < 6; ++i) error = std::fmax(error, std::fabs(r.x[i] - (1.0 + i)));
      std::printf("quadratic: %d cycles, %d evaluations, error %.2e\n", r.cycles, r.evaluations, error);
      check(r.converged && error < 1e-4, "ill-conditioned quadratic converges to the minimum");
   }

   {
      Vec x0(2); x0[0] = -1.2; x0[1] = 1.0;
      OptimizerResult r(Optimizer(tight).minimize(rosenbrock, x0));
      std::printf("rosenbrock: %d cycles, %d evaluations, x = (%.6f, %.6f)\n", r.cycles,
         r.evaluations, r.x[0], r.x[1]);
      check(r.converged && std::fabs(r.x[0] - 1.0) < 1e-4 && std::fabs(r.x[1] - 1.0) < 1e-4,
         "Rosenbrock function converges to (1,1)");
   }

   {
      // Pentagonal bipyramid, slightly distorted
      double ring = 1.0;
      Vec x0;
      for (int i = 0; i < 5; ++i) {
          double a(2.0 * M_PI * i / 5.0);
          x0.push_back(ring * std::cos(a) * 1.1); x0.push_back(ring * std::sin(a) * 0.95); x0.push_back(0.05 * i);
      }
      x0.push_back(0.0); x0.push_back(0.0); x0.push_back(0.9);
      x0.push_back(0.0); x0.push_back(0.0); x0.push_back(-0.9);
      OptimizerResult r(Optimizer(tight).minimize(lennardJones, x0));
      std::printf("LJ7: %d cycles, E = %.6f\n", r.cycles, r.energy);
      check(r.converged && std::fabs(r.energy + 16.505384) < 1e-4, "LJ7 cluster reaches the known minimum energy");
   }

   {
      // Default (xtb-like) criteria, and cancellation through the progress callback
      Vec x0(6, 0.0);
      OptimizerResult r(Optimizer().minimize(quadratic, x0));
      check(r.converged && r.gradientNorm < 1e-3, "default criteria are met");

      int calls(0);
      OptimizerResult c(Optimizer().minimize(quadratic, x0,
         [&calls](int, double, double, Vec const&) { return ++calls < 3; }));
      check(c.cancelled && !c.converged, "returning false from the progress callback cancels");
   }

   {
      // Failing evaluations are handled, and a bad start is reported
      int count(0);
      Vec x0(6, 0.0);
      OptimizerResult r(Optimizer().minimize([&count](Vec const& x, double& e, Vec& g) {
         if (++count == 3) return false;    // one failure in the middle of the run
         return quadratic(x, e, g);
      }, x0));
      check(r.converged, "a single failed evaluation is recovered from");

      OptimizerResult bad(Optimizer().minimize([](Vec const&, double&, Vec&) { return false; }, x0));
      check(!bad.converged && !bad.message.empty(), "failure at the start is reported");
   }

   {
      // Already at the minimum: converges immediately
      Vec x0; for (int i = 0; i < 6; ++i) x0.push_back(1.0 + i);
      OptimizerResult r(Optimizer().minimize(quadratic, x0));
      check(r.converged && r.cycles == 0, "a minimum is recognised without taking a step");
   }

   std::printf("%d failure(s)\n", s_failures);
   return s_failures == 0 ? 0 : 1;
}
