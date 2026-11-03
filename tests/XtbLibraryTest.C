// Copyright (C) 2026 Roberto Peverati.  Part of IQmol, free software under the GNU General Public
// License, version 3 or (at your option) any later version.  Developed with the assistance of
// Claude (Anthropic).
// Standalone test of libxtb (the C API of xtb): energies, gradients and solvation are
// compared with values from the xtb program for water.  It also runs an OpenMP region of
// our own first, because IQmol uses libomp while libxtb uses GNU libgomp and the two
// have to coexist in one process.
//
// Build:  make XtbLibraryTest      (needs -DXTB_ROOT=<xtb install dir> at configure time)

#include <xtb.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

static int s_failures(0);

static void check(bool ok, char const* what)
{
   std::printf("%s: %s\n", ok ? "PASS" : "FAIL", what);
   if (!ok) ++s_failures;
}

static double const AngstromToBohr(1.0 / 0.529177210903);

enum Method { GFN1, GFN2, GFNFF };

struct Result {
   bool   ok;
   double energy;
   double gsolv;
   std::vector<double> gradient;   // Hartree/Bohr
   std::string message;
};

// Water, atoms in Angstrom as in the xtb program reference calculations
static std::vector<double> water()
{
   double a[] = { 0.0, 0.0, 0.0,   0.0, 0.9, 0.6,   0.0, -0.9, 0.6 };
   return std::vector<double>(a, a + 9);
}

static Result run(std::vector<double> const& angstrom, Method method, double charge, int uhf,
   char const* solvent = 0, char const* model = "alpb")
{
   Result r;
   r.ok = false; r.energy = 0.0; r.gsolv = 0.0;
   int const natoms(3);
   int numbers[] = { 8, 1, 1 };

   std::vector<double> bohr(angstrom);
   for (size_t i = 0; i < bohr.size(); ++i) bohr[i] *= AngstromToBohr;

   xtb_TEnvironment env(xtb_newEnvironment());
   char const* outputFile("xtb_library_test_output.txt");
   if (solvent) {
      // xtb_getSolvationEnergy only works for CPCM-X, so the ALPB solvation free
      // energy has to be read from the program output
      xtb_setOutput(env, outputFile);
      xtb_setVerbosity(env, XTB_VERBOSITY_FULL);
   }else {
      xtb_setVerbosity(env, XTB_VERBOSITY_MUTED);
   }
   xtb_TMolecule mol(xtb_newMolecule(env, &natoms, numbers, bohr.data(), &charge, &uhf, 0, 0));
   xtb_TCalculator calc(xtb_newCalculator());
   xtb_TResults res(xtb_newResults());

   if (xtb_checkEnvironment(env)) { r.message = "setup failed"; }
   else {
      switch (method) {
         case GFN1:  xtb_loadGFN1xTB(env, mol, calc, 0); break;
         case GFN2:  xtb_loadGFN2xTB(env, mol, calc, 0); break;
         case GFNFF: xtb_loadGFNFF(env, mol, calc, 0);   break;
      }
      if (solvent) {
         // Without the model argument the library defaults to GBSA, not ALPB
         char name[64], mod[16];
         std::strncpy(name, solvent, 63); name[63] = 0;
         std::strncpy(mod, model, 15); mod[15] = 0;
         xtb_setSolvent(env, calc, name, 0, 0, 0, mod);
      }
      if (!xtb_checkEnvironment(env)) xtb_singlepoint(env, mol, calc, res);

      if (xtb_checkEnvironment(env)) {
         char buffer[512];
         int n(sizeof(buffer));
         xtb_getError(env, buffer, &n);
         r.message = buffer;
      }else {
         xtb_getEnergy(env, res, &r.energy);
         r.gradient.assign(9, 0.0);
         xtb_getGradient(env, res, r.gradient.data());
         if (solvent) {
            xtb_releaseOutput(env);
            std::ifstream in(outputFile);
            std::string line;
            while (std::getline(in, line)) {
               std::size_t k(line.find("-> Gsolv"));
               if (k != std::string::npos) r.gsolv = std::atof(line.c_str() + k + 8);
            }
            std::remove(outputFile);
         }
         r.ok = true;
      }
   }

   xtb_delResults(&res);
   xtb_delCalculator(&calc);
   xtb_delMolecule(&mol);
   xtb_delEnvironment(&env);
   return r;
}


int main()
{
   std::printf("xtb API version %d\n", xtb_getAPIVersion());

   // IQmol's own OpenMP runtime (libomp) is active before libxtb (libgomp) is used
#ifdef _OPENMP
   double sum(0.0);
   #pragma omp parallel for reduction(+:sum)
   for (int i = 0; i < 100000; ++i) sum += 1.0;
   check(sum == 100000.0, "OpenMP region of our own runs before libxtb is used");
   std::printf("OpenMP threads available: %d\n", omp_get_max_threads());
#else
   std::printf("(compiled without OpenMP)\n");
#endif

   double const tol(1e-6);   // Hartree

   Result gfn2(run(water(), GFN2, 0.0, 0));
   check(gfn2.ok, "GFN2-xTB single point runs");
   check(std::fabs(gfn2.energy + 5.049828064206) < tol, "GFN2-xTB energy matches the xtb program");

   Result gfn1(run(water(), GFN1, 0.0, 0));
   check(gfn1.ok && std::fabs(gfn1.energy + 5.746909443229) < tol, "GFN1-xTB energy matches the program");

   Result ff(run(water(), GFNFF, 0.0, 0));
   check(ff.ok && std::fabs(ff.energy + 0.309068176965) < tol, "GFN-FF energy matches the program");

   Result cation(run(water(), GFN2, 1.0, 1));
   check(cation.ok && std::fabs(cation.energy + 4.397622482702) < tol, "charged open-shell energy matches");

   Result alpb(run(water(), GFN2, 0.0, 0, "water"));
   check(alpb.ok && std::fabs(alpb.energy + 5.065761839366) < tol, "ALPB(water) energy matches");
   check(alpb.ok && std::fabs(alpb.gsolv + 0.022019795665) < tol, "ALPB(water) Gsolv matches");

   Result meoh(run(water(), GFN1, 0.0, 0, "methanol"));
   check(meoh.ok && std::fabs(meoh.energy + 5.758362579834) < tol, "GFN1 with ALPB(methanol) matches");

   // Gradient by central finite differences of the energy
   if (gfn2.ok) {
      double maxError(0.0);
      for (int i = 0; i < 9; ++i) {
          std::vector<double> plus(water()), minus(water());
          plus[i]  += 1e-3;    // Angstrom
          minus[i] -= 1e-3;
          Result p(run(plus, GFN2, 0.0, 0)), m(run(minus, GFN2, 0.0, 0));
          double numeric((p.energy - m.energy) / (2e-3 * AngstromToBohr));   // Eh/Bohr
          maxError = std::fmax(maxError, std::fabs(numeric - gfn2.gradient[i]));
      }
      std::printf("max |analytic - numeric gradient| = %.2e Eh/Bohr\n", maxError);
      check(maxError < 1e-4, "analytic gradient agrees with finite differences");
   }

   // The library does NOT validate the spin: 8 electrons with one 'unpaired' electron
   // silently gives the singlet energy, so IQmol has to check the parity itself.
   Result inconsistent(run(water(), GFN2, 0.0, 1));
   check(inconsistent.ok && std::fabs(inconsistent.energy - gfn2.energy) < tol,
      "inconsistent spin is NOT rejected by the library (documented, IQmol must check)");

   Result badSolvent(run(water(), GFN2, 0.0, 0, "notasolvent"));
   check(!badSolvent.ok, "unknown solvent is reported as an error");
   if (!badSolvent.ok) std::printf("  (message: %s)\n", badSolvent.message.c_str());

   // Our own OpenMP regions and xtb calculations (libgomp) alternating in one process
   bool stable(true);
   for (int i = 0; i < 25 && stable; ++i) {
#ifdef _OPENMP
      double t(0.0);
      #pragma omp parallel for reduction(+:t)
      for (int k = 0; k < 10000; ++k) t += 1.0;
      if (t != 10000.0) stable = false;
#endif
      Result again(run(water(), GFN2, 0.0, 0));
      if (!again.ok || std::fabs(again.energy - gfn2.energy) > 1e-9) stable = false;
   }
   check(stable, "25 alternating OpenMP regions and xtb calculations are stable and reproducible");

   std::printf("%d failure(s)\n", s_failures);
   return s_failures == 0 ? 0 : 1;
}
