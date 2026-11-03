# Copyright (C) 2026 Roberto Peverati
# Developed with the assistance of Claude (Anthropic).
#
# This file is part of IQMolX, a fork of IQmol (http://iqmol.org). It is free
# software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version. It is
# distributed WITHOUT ANY WARRANTY; see <http://www.gnu.org/licenses/>.
"""Example plugin for IQMolX: Lennard-Jones noble gas clusters.

An engine is a subclass of iqmolx_engines.base.Engine.  The class attributes describe it (this
is everything IQMolX shows in its window), load() prepares it and compute() does the work.
"""

import numpy as np

from iqmolx_engines.base import Engine, option

# Rounded textbook parameters: sigma in Angstrom, epsilon in eV
PARAMETERS = {
    2:  (2.56, 0.00088),     # He
    10: (2.75, 0.00307),     # Ne
    18: (3.40, 0.01032),     # Ar
    36: (3.60, 0.01400),     # Kr
    54: (4.10, 0.01980),     # Xe
}


class LennardJonesEngine(Engine):
    name = "lj"                                    # what load requests call it
    title = "Lennard-Jones (example plugin)"       # what the window shows
    description = "Noble gas clusters with a 12-6 pair potential. An example, not a research method."
    license = "GPL-3.0-or-later"
    elements = tuple(sorted(PARAMETERS))           # He, Ne, Ar, Kr, Xe
    charge = False                                 # neutral atoms only
    multiplicity = False                           # closed shell only
    properties = ("energy", "forces")              # energy and forces are mandatory
    options = (
        # scope "load": applied when the model loads (a change loads it again)
        option("epsilon_scale", "Well depth scale", "float", 1.0, min=0.1, max=10.0,
               description="Multiplies every epsilon."),
    )

    def load(self, options, ctx):
        # Nothing to load here.  A real engine loads its model (and may download weights, with
        # ctx.progress(fraction, message) reporting how far it is) and returns a dict of
        # information for the status line.
        self.scale = options["epsilon_scale"]
        return {"model": "Lennard-Jones"}

    def compute(self, numbers, coordinates, charge, multiplicity, want, options, ctx):
        n = len(numbers)
        sigma = np.array([PARAMETERS[int(z)][0] for z in numbers])
        eps = np.array([PARAMETERS[int(z)][1] for z in numbers]) * self.scale

        energy = 0.0
        forces = np.zeros((n, 3))
        for i in range(n - 1):
            ctx.check_cancel()                     # a long computation lets the user cancel
            d = coordinates[i] - coordinates[i + 1:]
            r = np.linalg.norm(d, axis=1)
            s = 0.5 * (sigma[i] + sigma[i + 1:])             # Lorentz mixing
            e = np.sqrt(eps[i] * eps[i + 1:])                # Berthelot mixing
            x6 = (s / r) ** 6
            energy += float(np.sum(4.0 * e * (x6 * x6 - x6)))
            # -dE/dr along the pair, a repulsion when positive
            f = 24.0 * e * (2.0 * x6 * x6 - x6) / r
            pair = (f / r)[:, None] * d
            forces[i] += pair.sum(axis=0)
            forces[i + 1:] -= pair

        # Return what was asked for (extra keys are ignored); units are eV and eV/Angstrom
        return {"energy": energy, "forces": forces}
