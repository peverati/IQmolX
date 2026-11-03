# Copyright (C) 2026 Roberto Peverati
# Developed with the assistance of Claude (Anthropic).
#
# This file is part of IQMolX, a fork of IQmol (http://iqmol.org). It is free
# software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version. It is
# distributed WITHOUT ANY WARRANTY; see <http://www.gnu.org/licenses/>.

"""Demo engine: a pairwise Morse potential. It is a toy, not a real method; it exists so
the protocol, the IQMolX window and the tests have an engine with no dependencies except
numpy. It has analytic forces and Hessian, a dipole from equally shared charge, and an
optional artificial delay that polls for cancellation."""

import time

import numpy as np

from ..base import Engine, option

# Covalent radii in Angstrom (Cordero et al., Dalton Trans. 2008, 2832), H to Ar.
RADII = {1: 0.31, 2: 0.28, 3: 1.28, 4: 0.96, 5: 0.84, 6: 0.76, 7: 0.71, 8: 0.66, 9: 0.57,
         10: 0.58, 11: 1.66, 12: 1.41, 13: 1.21, 14: 1.11, 15: 1.07, 16: 1.05, 17: 1.02,
         18: 1.06}


def morse(numbers, xyz, depth, alpha, want_hessian=False):
    """Energy, forces and optionally the Hessian of V = sum_pairs D[(1-e)^2 - 1] with
    e = exp(-alpha (r - r0)) and r0 the sum of covalent radii."""
    n = len(numbers)
    energy = 0.0
    forces = np.zeros((n, 3))
    hess = np.zeros((3 * n, 3 * n)) if want_hessian else None
    eye = np.eye(3)
    for i in range(n):
        for j in range(i + 1, n):
            d = xyz[i] - xyz[j]
            r = float(np.linalg.norm(d))
            r0 = RADII[int(numbers[i])] + RADII[int(numbers[j])]
            e = np.exp(-alpha * (r - r0))
            energy += depth * ((1.0 - e) ** 2 - 1.0)
            dv = 2.0 * depth * alpha * e * (1.0 - e)
            u = d / r
            forces[i] -= dv * u
            forces[j] += dv * u
            if want_hessian:
                d2v = 2.0 * depth * alpha ** 2 * e * (2.0 * e - 1.0)
                uu = np.outer(u, u)
                block = d2v * uu + (dv / r) * (eye - uu)
                a, b = slice(3 * i, 3 * i + 3), slice(3 * j, 3 * j + 3)
                hess[a, a] += block
                hess[b, b] += block
                hess[a, b] -= block
                hess[b, a] -= block
    return energy, forces, hess


class DemoEngine(Engine):
    name = "demo"
    title = "Demo (Morse toy potential)"
    description = ("Pairwise Morse potential for tests and demonstrations. "
                   "Not a real method: energies have no chemical meaning.")
    license = "GPL-3.0-or-later"
    elements = tuple(sorted(RADII))
    charge = True
    multiplicity = True
    properties = ("energy", "forces", "hessian", "dipole", "charges")
    options = (
        option("well_depth", "Well depth (eV)", "float", 4.0, min=0.1, max=20.0,
               description="Depth D of every pair potential."),
        option("alpha", "Width parameter (1/Angstrom)", "float", 1.8, min=0.2, max=5.0,
               description="Steepness of the Morse potential."),
        option("delay_ms", "Artificial delay (ms)", "int", 0, scope="compute", min=0,
               max=60000, description="Sleep before answering; polls for cancellation."),
    )

    def load(self, options, ctx):
        self.depth = options["well_depth"]
        self.alpha = options["alpha"]
        return {"device": "cpu", "model": "morse"}

    def compute(self, numbers, coordinates, charge, multiplicity, want, options, ctx):
        deadline = time.monotonic() + options["delay_ms"] / 1000.0
        while time.monotonic() < deadline:
            ctx.check_cancel()
            time.sleep(0.005)
        ctx.check_cancel()
        energy, forces, hess = morse(numbers, coordinates, self.depth, self.alpha,
                                     "hessian" in want)
        result = {"energy": energy, "forces": forces}
        if hess is not None:
            result["hessian"] = hess
        q = np.full(len(numbers), charge / len(numbers))
        if "charges" in want:
            result["charges"] = q
        if "dipole" in want:
            result["dipole"] = q @ coordinates
        return result
