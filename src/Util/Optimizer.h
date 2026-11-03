#pragma once
/*******************************************************************************

  This file is part of IQmol, a free molecular visualization program. See
  <http://iqmol.org> for more details.

  IQmol is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

********************************************************************************/

#include <cmath>
#include <functional>
#include <string>
#include <vector>

/*!
 *  \file Optimizer.h
 *
 *  \brief Limited-memory BFGS minimizer for functions of many variables, used to
 *  optimize molecular geometries with backends that only provide energies and
 *  gradients (such as libxtb).  It has no dependencies on Qt or on the backend.
 *
 *  The convergence criteria follow xtb's 'normal' level: the energy change of
 *  the last step is below 5e-6 Eh and the gradient norm below 1e-3 Eh/Bohr.
 */

namespace IQmol {
namespace Calc {

   struct OptimizerSettings {
      int    maxCycles;     // accepted steps
      double energyTol;     // Hartree
      double gradientTol;   // norm of the gradient, Hartree/Bohr
      double maxStep;       // largest change of any one coordinate in a step, Bohr
      int    history;       // number of correction pairs kept

      OptimizerSettings() : maxCycles(500), energyTol(5e-6), gradientTol(1e-3),
         maxStep(0.3), history(10) { }
   };

   struct OptimizerResult {
      bool   converged;
      bool   cancelled;
      int    cycles;        // number of accepted steps
      int    evaluations;   // number of energy/gradient evaluations
      double energy;
      double gradientNorm;
      std::vector<double> x;
      std::string message;

      OptimizerResult() : converged(false), cancelled(false), cycles(0), evaluations(0),
         energy(0.0), gradientNorm(0.0) { }
   };


   class Optimizer {

      public:
         /// Computes the energy and gradient at x.  Returns false if that fails.
         typedef std::function<bool(std::vector<double> const& x, double& energy,
            std::vector<double>& gradient)> Evaluate;

         /// Called before each step with the current point.  Return false to stop.
         typedef std::function<bool(int cycle, double energy, double gradientNorm,
            std::vector<double> const& x)> Progress;

         explicit Optimizer(OptimizerSettings const& settings = OptimizerSettings())
          : m_settings(settings) { }

         OptimizerResult minimize(Evaluate const& evaluate, std::vector<double> const& x0,
            Progress const& progress = Progress()) const
         {
            OptimizerResult result;
            std::vector<double> x(x0), g, xNew, gNew;
            double energy(0.0), energyNew(0.0);
            size_t const n(x.size());

            ++result.evaluations;
            if (!evaluate(x, energy, g)) {
               result.message = "Energy evaluation failed at the starting geometry";
               result.x = x;
               return result;
            }

            std::vector<std::vector<double> > S, Y;   // correction pairs, oldest first
            std::vector<double> rho;
            double gamma(1.0);                        // initial inverse Hessian scale, Bohr^2/Eh
            double lastChange(1.0);                   // energy change of the last step

            for (int cycle = 0; ; ++cycle) {
                double gnorm(norm(g));
                result.energy = energy;
                result.gradientNorm = gnorm;
                result.cycles = cycle;
                result.x = x;

                if (progress && !progress(cycle, energy, gnorm, x)) {
                   result.cancelled = true;
                   result.message = "Optimization cancelled";
                   return result;
                }

                if (gnorm < m_settings.gradientTol &&
                    (cycle == 0 || std::fabs(lastChange) < m_settings.energyTol)) {
                   result.converged = true;
                   return result;
                }

                if (cycle >= m_settings.maxCycles) {
                   result.message = "Maximum number of optimization cycles reached";
                   return result;
                }

                // Search direction from the two-loop recursion
                std::vector<double> d(direction(g, S, Y, rho, gamma));
                if (dot(d, g) >= 0.0) {           // not a descent direction: start again
                   S.clear(); Y.clear(); rho.clear();
                   d = direction(g, S, Y, rho, gamma);
                }

                // Limit the largest coordinate change
                double dmax(0.0);
                for (size_t i = 0; i < n; ++i) dmax = std::fmax(dmax, std::fabs(d[i]));
                if (dmax > m_settings.maxStep) {
                   double scale(m_settings.maxStep / dmax);
                   for (size_t i = 0; i < n; ++i) d[i] *= scale;
                }

                // Backtracking line search on the energy (Armijo condition)
                double slope(dot(d, g));
                double t(1.0);
                bool accepted(false);
                for (int trial = 0; trial < 10; ++trial) {
                    xNew = x;
                    for (size_t i = 0; i < n; ++i) xNew[i] += t * d[i];
                    ++result.evaluations;
                    bool ok(evaluate(xNew, energyNew, gNew));
                    if (ok && energyNew <= energy + 1e-4 * t * slope) { accepted = true; break; }
                    t *= 0.5;
                }

                if (!accepted) {
                   if (!S.empty()) {            // the curvature model is bad: reset it
                      S.clear(); Y.clear(); rho.clear();
                      gamma = 1.0;
                      continue;
                   }
                   result.message = "Line search failed to lower the energy";
                   return result;
                }

                // Update the curvature information
                std::vector<double> s(n), y(n);
                for (size_t i = 0; i < n; ++i) { s[i] = xNew[i] - x[i]; y[i] = gNew[i] - g[i]; }
                double sy(dot(s, y));
                if (sy > 1e-10) {
                   if (static_cast<int>(S.size()) >= m_settings.history) {
                      S.erase(S.begin()); Y.erase(Y.begin()); rho.erase(rho.begin());
                   }
                   S.push_back(s); Y.push_back(y); rho.push_back(1.0 / sy);
                   gamma = sy / dot(y, y);
                }

                lastChange = energyNew - energy;
                x = xNew; g = gNew; energy = energyNew;
            }
         }

      private:
         static double dot(std::vector<double> const& a, std::vector<double> const& b)
         {
            double sum(0.0);
            for (size_t i = 0; i < a.size(); ++i) sum += a[i] * b[i];
            return sum;
         }

         static double norm(std::vector<double> const& a) { return std::sqrt(dot(a, a)); }

         // Standard L-BFGS two-loop recursion: returns -H*g
         static std::vector<double> direction(std::vector<double> const& g,
            std::vector<std::vector<double> > const& S, std::vector<std::vector<double> > const& Y,
            std::vector<double> const& rho, double gamma)
         {
            size_t const m(S.size()), n(g.size());
            std::vector<double> q(g), alpha(m);
            for (size_t k = m; k-- > 0; ) {
                alpha[k] = rho[k] * dot(S[k], q);
                for (size_t i = 0; i < n; ++i) q[i] -= alpha[k] * Y[k][i];
            }
            for (size_t i = 0; i < n; ++i) q[i] *= gamma;
            for (size_t k = 0; k < m; ++k) {
                double beta(rho[k] * dot(Y[k], q));
                for (size_t i = 0; i < n; ++i) q[i] += S[k][i] * (alpha[k] - beta);
            }
            for (size_t i = 0; i < n; ++i) q[i] = -q[i];
            return q;
         }

         OptimizerSettings m_settings;
   };

} } // end namespace IQmol::Calc
