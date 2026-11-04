# Copyright (C) 2026 Roberto Peverati
# Developed with the assistance of Claude (Anthropic).
#
# This file is part of IQMolX, a fork of IQmol (http://iqmol.org). It is free
# software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version. It is
# distributed WITHOUT ANY WARRANTY; see <http://www.gnu.org/licenses/>.

"""Example plugin for IQMolX: the effective medium theory (EMT) calculator of ASE.

This is what turning an ASE calculator into an engine of its own looks like.  The generic "ASE
calculator" engine of IQMolX can already run EMT (type ase.calculators.emt:EMT in its window), but
a plugin of its own gives the calculator a proper name and description, only the elements it can
do, the right license, and settings with a meaning in the window instead of JSON text.
"""

import importlib.util

import numpy as np

from iqmolx_engines.base import Engine, EngineError, option


class EmtEngine(Engine):
    name = "emt"                                    # what load requests call it
    title = "EMT metals (ASE, example plugin)"      # what the window shows
    description = ("Effective medium theory (Jacobsen, Norskov, Puska 1987) as implemented in ASE: a "
                   "very fast empirical potential for the metals Al, Ni, Cu, Pd, Ag, Pt and Au, with "
                   "H, C, N and O for small adsorbates.  Rough, meant for tests and for learning.")
    license = "LGPL-2.1-or-later (ASE)"
    citation = "Jacobsen, Stoltze, Norskov, Surf. Sci. 366, 394 (1996)"
    elements = (1, 6, 7, 8, 13, 28, 29, 46, 47, 78, 79)       # H C N O Al Ni Cu Pd Ag Pt Au
    charge = False                                  # neutral only
    multiplicity = False                            # no spin
    properties = ("energy", "forces")               # energy and forces are mandatory
    requirements = ("ase",)                         # what the Install button of the window installs
    options = (
        # scope "load" (the default): applied when the model loads
        option("asap_cutoff", "ASAP cutoff", "bool", False,
               description="Use the cutoff of the ASAP program instead of the ASE one; the energies "
                           "differ by about 1e-6 eV."),
        # scope "compute": sent with every calculation, so changing it needs no new load
        option("min_distance", "Closest atoms (Angstrom)", "float", 0.5, scope="compute", min=0.0, max=2.0,
               description="EMT is meaningless when atoms almost touch (two copper atoms at 0.05 "
                           "Angstrom have 1264 eV).  A structure with atoms closer than this is refused; "
                           "0 switches the check off."),
    )

    @classmethod
    def unavailable_reason(cls):
        # The engine is listed even when ASE is missing: the window then says what is missing and
        # offers to install it.  Do not import ASE at the top of the module for this reason.
        if importlib.util.find_spec("ase") is None:
            return "ASE is not installed (pip install ase)"
        return ""

    def load(self, options, ctx):
        from ase.calculators.emt import EMT
        self.calc = EMT(asap_cutoff=options["asap_cutoff"])
        return {"calculator": "ase.calculators.emt.EMT"}       # shown in the status line

    def compute(self, numbers, coordinates, charge, multiplicity, want, options, ctx):
        from ase import Atoms
        atoms = Atoms(numbers=numbers, positions=coordinates)      # Angstrom, no periodic cell
        if len(numbers) > 1 and options["min_distance"] > 0.0:
            distances = atoms.get_all_distances()
            np.fill_diagonal(distances, np.inf)
            if distances.min() < options["min_distance"]:
                raise EngineError("engine_error", "Two atoms are only %.2f Angstrom apart; EMT is not "
                                  "meaningful for that" % distances.min(),
                                  detail="Move the atoms apart, or lower the closest distance allowed "
                                         "in the settings.")
        atoms.calc = self.calc
        # Units are eV and eV/Angstrom, the units of ASE and of the protocol: no conversion
        return {"energy": float(atoms.get_potential_energy()),
                "forces": np.asarray(atoms.get_forces())}
