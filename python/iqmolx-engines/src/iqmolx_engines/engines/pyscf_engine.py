# Copyright (C) 2026 Roberto Peverati
# Developed with the assistance of Claude (Anthropic).
#
# This file is part of IQMolX, a fork of IQmol (http://iqmol.org). It is free
# software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version. It is
# distributed WITHOUT ANY WARRANTY; see <http://www.gnu.org/licenses/>.
"""PySCF engine: Hartree-Fock and density functional theory with PySCF (Apache-2.0), with
analytic energies, forces and Hessians, a dipole and Mulliken charges.  Any charge and spin
multiplicity (restricted for closed shells, unrestricted for open shells), any Gaussian basis set that
PySCF knows.  This is an ab initio method: much slower than a neural network potential, practical for
small molecules.  PySCF is imported only when a calculation runs."""

import importlib.util
import math
import sys

import numpy as np

from ..base import Engine, EngineError, option

HARTREE_TO_EV = 27.211386245988
BOHR_IN_ANGSTROM = 0.529177210903
DEBYE_TO_E_ANGSTROM = 0.20819434          # 1 debye in e*Angstrom

# The menu: functionals that were run through energy, gradient and Hessian with PySCF 2.14, and
# "other", which opens a field for the name of any functional PySCF (libxc) knows.  wB97X-D is not
# implemented in PySCF and the -D3 variants need the extra package pyscf-dispersion.
METHODS = ["hf", "b3lyp", "pbe0", "pbe", "bp86", "tpss", "m06-2x", "wb97x", "r2scan", "wb97m-v", "other"]

# For meta-GGA functionals the analytic Hessian of PySCF 2.14 is not accurate (errors of 0.6 to 6
# eV/Angstrom^2 for water, against 0.05 for the others), so the Hessian is built from differences of
# the analytic gradients instead: slower (6N gradients), but correct.  PySCF itself says which
# functionals are meta-GGAs (libxc.is_meta_gga), so this also holds for a functional typed in "other".
STEP = 0.005                                  # Angstrom, for those finite differences

class PySCFEngine(Engine):
    name = "pyscf"
    title = "PySCF (Hartree-Fock and DFT)"
    description = ("Hartree-Fock and density functional theory with PySCF: analytic energies, forces and "
                   "Hessians for any charge and spin multiplicity, with the Gaussian basis set of your "
                   "choice. An ab initio method, much slower than a neural network potential: practical "
                   "for small molecules. For meta-GGA functionals (tpss, m06-2x, r2scan, wb97m-v and others) the "
                   "Hessian is built from differences of forces, which takes about 6N times as long.")
    license = "Apache-2.0 (PySCF)"
    citation = ("Sun et al., J. Chem. Phys. 153, 024109 (2020); "
                "Sun, WIREs Comput. Mol. Sci. 8, e1340 (2018)")
    elements = None
    charge = True
    multiplicity = True
    properties = ("energy", "forces", "hessian", "dipole", "charges")
    requirements = ("pyscf",)
    options = (
        option("method", "Method", "choice", "b3lyp", choices=METHODS,
               description="Hartree-Fock (hf) or a density functional. wb97m-v is accurate and "
                           "slower (it has a nonlocal VV10 part); other lets you type the name of "
                           "any functional that PySCF knows."),
        option("method_other", "Functional", "string", "", required=True, placeholder="cam-b3lyp",
               show_if={"method": "other"},
               description="The name of a functional of PySCF (libxc), for example cam-b3lyp, wb97x-v, "
                           "scan or pbe0, or hf. Names that end in -d3bj need the package "
                           "pyscf-dispersion."),
        option("basis", "Basis set", "string", "def2-svp", required=True, placeholder="def2-svp",
               description="Any basis set of PySCF, for example sto-3g, 6-31g*, def2-svp, def2-tzvp, "
                           "cc-pvdz, aug-cc-pvdz."),
        option("density_fit", "Density fitting", "bool", True,
               description="Speeds up the calculation, with a small error (about 1e-4 Hartree for water)."),
        option("grid_level", "Integration grid (DFT)", "int", 3, min=0, max=9,
               description="Level of the numerical integration grid of the density functionals; "
                           "3 is the PySCF default, higher is finer and slower."),
        option("nlc_grid_level", "Grid for the VV10 part", "int", 1, min=0, max=9,
               show_if={"method": "wb97m-v"},
               description="Level of the grid of the nonlocal (VV10) term of functionals such as "
                           "wb97m-v. Level 1 gives energies within 1e-6 Hartree of level 3 (the "
                           "PySCF default) and is about 8 times faster; 0 is faster still, 2 is finer."),
        option("max_cycle", "Maximum SCF cycles", "int", 100, min=10, max=1000),
    )

    @classmethod
    def unavailable_reason(cls):
        if importlib.util.find_spec("pyscf") is None:
            if sys.platform == "win32":
                return ("pyscf not installed. The PySCF project does not support native Windows "
                        "(it advises the Windows Subsystem for Linux); 'pip install pyscf' may "
                        "work in a Python of your own")
            return "pyscf not installed (pip install pyscf)"
        return ""

    def load(self, options, ctx):
        import pyscf
        from pyscf.dft import libxc

        self.settings = dict(options)
        method = options["method_other"] if options["method"] == "other" else options["method"]
        self.xc = method.strip().lower()
        if self.xc != "hf":
            try:
                libxc.parse_xc(self.xc)
            except Exception as exc:
                raise EngineError("bad_option", "PySCF does not know the functional '%s'" % method.strip(),
                                  detail="%s\nUse a name of libxc, for example cam-b3lyp, wb97x-v, scan." % exc) from None
        self.meta_gga = self.xc != "hf" and bool(libxc.is_meta_gga(self.xc))
        self.nonlocal_term = self.xc != "hf" and bool(libxc.nlc_coeff(self.xc))
        self.last_key = None          # the system of the last calculation, and its density
        self.last_density = None
        info = {"method": self.xc, "basis": options["basis"], "density_fit": options["density_fit"],
                "pyscf": pyscf.__version__}
        if self.meta_gga:
            info["hessian"] = "from differences of forces (meta-GGA)"
        return info

    # -- helpers -------------------------------------------------------------------------------

    @staticmethod
    def check_spin(numbers, charge, multiplicity):
        electrons = int(sum(int(z) for z in numbers)) - charge
        unpaired = multiplicity - 1
        if electrons < 0 or unpaired > electrons or (electrons - unpaired) % 2:
            raise EngineError("unsupported_multiplicity",
                              "A charge of %d and a spin multiplicity of %d are not possible for %d electrons"
                              % (charge, multiplicity, electrons),
                              detail="An even number of electrons needs an odd multiplicity (singlet, "
                                     "triplet, ...) and an odd number of electrons an even one (doublet, ...).")

    def build_molecule(self, numbers, coordinates, charge, multiplicity):
        from pyscf import gto
        mol = gto.Mole()
        mol.atom = [[int(z), tuple(float(c) for c in xyz)] for z, xyz in zip(numbers, coordinates)]
        mol.unit = "Angstrom"
        mol.basis = self.settings["basis"].strip()
        mol.charge = int(charge)
        mol.spin = int(multiplicity) - 1
        mol.symmetry = False
        mol.verbose = 0
        mol.output = None
        try:
            mol.build()
        except Exception as exc:
            raise EngineError("engine_error", "PySCF cannot set up the molecule with the basis '%s': %s"
                              % (self.settings["basis"], exc),
                              detail="Check the spelling of the basis set and that it has every element "
                                     "of the molecule.") from None
        return mol

    def build_method(self, mol, ctx):
        from pyscf import dft, scf
        s = self.settings
        open_shell = mol.spin != 0
        if self.xc == "hf":
            mf = scf.UHF(mol) if open_shell else scf.RHF(mol)
        else:
            mf = dft.UKS(mol) if open_shell else dft.RKS(mol)
            mf.xc = self.xc
            mf.grids.level = s["grid_level"]
            if self.nonlocal_term:
                mf.nlcgrids.level = s["nlc_grid_level"]
        if s["density_fit"]:
            mf = mf.density_fit()
        mf.max_cycle = s["max_cycle"]
        mf.conv_tol = 1e-9
        mf.callback = lambda envs: ctx.check_cancel()          # the SCF can be cancelled between cycles
        return mf

    def run_scf(self, mf, key):
        guess = self.last_density if key == self.last_key else None
        mf.kernel(dm0=guess)
        if not mf.converged:                                    # one more try from the default guess
            mf.level_shift = 0.3
            mf.kernel()
        if not mf.converged:
            raise EngineError("engine_error", "The SCF did not converge",
                              detail="Try another geometry, a smaller basis set, or no density fitting.")
        self.last_key, self.last_density = key, mf.make_rdm1()

    # -- the calculation -------------------------------------------------------------------------

    def scf(self, numbers, coordinates, charge, multiplicity, ctx):
        """Converged SCF of one geometry."""
        mol = self.build_molecule(numbers, coordinates, charge, multiplicity)
        mf = self.build_method(mol, ctx)
        self.run_scf(mf, (tuple(int(z) for z in numbers), int(charge), int(multiplicity)))
        return mf

    def forces_of(self, mf):
        """Forces in eV/Angstrom."""
        method = mf.nuc_grad_method()
        if self.xc != "hf":
            method.grid_response = True            # the integration grid moves with the atoms: more accurate
        return -np.asarray(method.kernel()) * HARTREE_TO_EV / BOHR_IN_ANGSTROM       # kernel: Hartree/Bohr

    def analytic_hessian(self, mf, n):
        h = np.asarray(mf.Hessian().kernel())                              # (N, N, 3, 3), Hartree/Bohr^2
        h = h.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n)
        return h * HARTREE_TO_EV / BOHR_IN_ANGSTROM ** 2

    def numerical_hessian(self, numbers, coordinates, charge, multiplicity, ctx):
        """eV/Angstrom^2 from central differences of the analytic forces."""
        n = len(numbers)
        x0 = np.asarray(coordinates, dtype=float)
        h = np.zeros((3 * n, 3 * n))
        for i in range(3 * n):
            ctx.check_cancel()
            forces = []
            for sign in (+1, -1):
                x = x0.copy()
                x[i // 3, i % 3] += sign * STEP
                forces.append(self.forces_of(self.scf(numbers, x, charge, multiplicity, ctx)).ravel())
            h[i] = -(forces[0] - forces[1]) / (2 * STEP)
        return h

    def compute(self, numbers, coordinates, charge, multiplicity, want, options, ctx):
        ctx.check_cancel()
        n = len(numbers)
        self.check_spin(numbers, charge, multiplicity)
        mf = self.scf(numbers, coordinates, charge, multiplicity, ctx)

        result = {"energy": float(mf.e_tot) * HARTREE_TO_EV}
        if "forces" in want:
            ctx.check_cancel()
            result["forces"] = self.forces_of(mf)
        if "dipole" in want:
            result["dipole"] = np.asarray(mf.dip_moment(verbose=0)) * DEBYE_TO_E_ANGSTROM
        if "charges" in want:
            result["charges"] = np.asarray(mf.mulliken_pop(verbose=0)[1])
        if "hessian" in want:
            ctx.check_cancel()
            if self.meta_gga:
                h = self.numerical_hessian(numbers, coordinates, charge, multiplicity, ctx)
            else:
                h = self.analytic_hessian(mf, n)
            result["hessian"] = 0.5 * (h + h.T)    # symmetric up to numerical noise (grid, differences)
        return result
