# Copyright (C) 2026 Roberto Peverati
# Developed with the assistance of Claude (Anthropic).
#
# This file is part of IQMolX, a fork of IQmol (http://iqmol.org). It is free
# software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version. It is
# distributed WITHOUT ANY WARRANTY; see <http://www.gnu.org/licenses/>.
"""Egret engines: the Egret family of neural network potentials of Rowan Scientific, built on the
MACE architecture (https://github.com/rowansci/egret-public).  The code and the model weights are
MIT licensed (README and preprint of the repository), and so is mace-torch.  The models are for
neutral molecules of H, C, N, O, F, P, S, Cl, Br and I (Egret-1e also Si) in the gas phase; they
have no input for charge or spin.

The weights are single files in the Egret repository.  They are downloaded the first time a model is
used, from a pinned commit, checked against the hash that GitHub lists for the file, and kept in the
cache folder of the user (see download.py).  The models are run with the calculator of mace-torch,
which is imported only when a model is loaded."""

import importlib.util
import os
import sys

import numpy as np

from .. import download
from ..base import Engine, EngineError, option

# The commit of rowansci/egret-public that the weights come from.  A file is only taken from this
# commit, so what is downloaded cannot change under us.
COMMIT = "ed85f3c7a792b5426dd73925e7f346590b60b2ad"
URL = "https://raw.githubusercontent.com/rowansci/egret-public/%s/compiled_models/%s"

# H, C, N, O, F, P, S, Cl, Br, I
ELEMENTS = (1, 6, 7, 8, 9, 15, 16, 17, 35, 53)
ELEMENTS_1E = (1, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53)     # and Si

PAPER = ("Egret-1: Pretrained Neural Network Potentials for Efficient and Accurate Bioorganic "
         "Simulation, Rowan Scientific (preprint, 2025), "
         "rowansci.com/publications/egret-1-pretrained-neural-network-potentials")
MACE = ("Batatia et al., MACE: Higher order equivariant message passing neural networks for fast "
        "and accurate force fields, NeurIPS 2022")


class EgretBase(Engine):
    """What the Egret models share.  A model sets `name`, `title`, `description`, `elements`,
    and `model_file` with the size and git blob hash of that file in the pinned commit."""

    model_file = ""
    model_size = 0
    blob_sha1 = ""
    license = "MIT (the Egret models, Rowan Scientific; mace-torch is MIT)"
    citation = PAPER + "; " + MACE
    elements = ELEMENTS
    charge = False
    multiplicity = False
    # No Hessian: the models are TorchScript exports, whose Hessian routine is a Python function that
    # cannot run in the export ("annotated to be ignored").  IQMolX builds frequencies from forces.
    properties = ("energy", "forces")
    requires_download = True
    requirements = ("mace-torch>=0.3.14",)
    options = (
        option("precision", "Precision", "choice", "float64", choices=["float64", "float32"],
               description="float64 gives smoother forces for optimizations and frequencies."),
    )

    @classmethod
    def unavailable_reason(cls):
        if sys.version_info < (3, 11):
            return "mace-torch needs Python 3.11 or newer"
        if importlib.util.find_spec("mace") is None:
            return "mace-torch not installed (pip install mace-torch)"
        return ""

    @classmethod
    def weights_path(cls):
        return os.path.join(download.cache_dir(), "egret", COMMIT[:10], cls.model_file)

    def load(self, options, ctx):
        from mace.calculators import MACECalculator

        path = download.fetch(URL % (COMMIT, self.model_file), self.weights_path(), ctx,
                              label="the %s weights" % self.title, size=self.model_size,
                              blob_sha1=self.blob_sha1)
        ctx.check_cancel()
        ctx.progress(None, "Loading %s" % self.title)
        try:
            self.calc = MACECalculator(model_paths=path, device="cpu",
                                       default_dtype=options["precision"])
        except Exception as exc:
            raise EngineError("engine_error", "Cannot load %s: %s" % (self.title, exc),
                              detail="The model file is %s. Egret needs mace-torch 0.3.14 or "
                                     "newer; try: uv pip install -U mace-torch" % path) from None
        return {"device": "cpu", "model": self.title, "weights": os.path.basename(path)}

    def compute(self, numbers, coordinates, charge, multiplicity, want, options, ctx):
        from ase import Atoms

        ctx.check_cancel()
        atoms = Atoms(numbers=[int(z) for z in numbers], positions=np.asarray(coordinates, dtype=float))
        atoms.calc = self.calc
        return {"energy": float(atoms.get_potential_energy()),
                "forces": np.asarray(atoms.get_forces(), dtype=float)}


class Egret1Engine(EgretBase):
    name = "egret1"
    title = "Egret-1 (Rowan, MACE)"
    description = ("Egret-1: a general-purpose neural network potential optimized for bioorganic "
                   "molecules, built on the MACE architecture. Neutral closed-shell molecules of H, C, "
                   "N, O, F, P, S, Cl, Br and I in the gas phase. Downloads about 23 MB the first time.")
    model_file, model_size = "EGRET_1.model", 23411282
    blob_sha1 = "17e750f256622a4e8ebe6f29109793f0c0e1d431"


class Egret1eEngine(EgretBase):
    name = "egret1e"
    title = "Egret-1e (Rowan, MACE)"
    description = ("Egret-1e: Egret-1 with main-group chemistry data added; excels at thermochemistry. "
                   "Neutral closed-shell molecules of H, C, N, O, F, Si, P, S, Cl, Br and I in the gas "
                   "phase. Downloads about 24 MB the first time.")
    elements = ELEMENTS_1E
    model_file, model_size = "EGRET_1E.model", 24305728
    blob_sha1 = "c8fb2653eaa59c7d7bfcd77a89c622c1cffdbe84"


class Egret1tEngine(EgretBase):
    name = "egret1t"
    title = "Egret-1t (Rowan, MACE)"
    description = ("Egret-1t: Egret-1 with transition states added; meant for chemical reactivity. "
                   "Neutral closed-shell molecules of H, C, N, O, F, P, S, Cl, Br and I in the gas "
                   "phase. Downloads about 23 MB the first time.")
    model_file, model_size = "EGRET_1T.model", 23424084
    blob_sha1 = "7dc6ac390289d6118a96a2d9df219f31aff4401f"


class Egret1mEngine(EgretBase):
    name = "egret1m"
    title = "Egret-1M (Rowan, MACE)"
    description = ("Egret-1M: a smaller Egret-1 with low-resolution equivariance (128 spherical "
                   "channels); fast, with a smaller memory footprint. Neutral closed-shell molecules "
                   "of H, C, N, O, F, P, S, Cl, Br and I in the gas phase. Downloads about 9 MB the "
                   "first time.")
    model_file, model_size = "EGRET_1M.model", 9435130
    blob_sha1 = "9d067a321a4dd0dd638ec9ee79f1e64902e8a89b"


class Egret1sEngine(EgretBase):
    name = "egret1s"
    title = "Egret-1S (Rowan, MACE)"
    description = ("Egret-1S: a smaller Egret-1 that is invariant only (92 spherical channels); the "
                   "fastest Egret model with the smallest memory footprint. Neutral closed-shell "
                   "molecules of H, C, N, O, F, P, S, Cl, Br and I in the gas phase. Downloads about "
                   "4 MB the first time.")
    model_file, model_size = "EGRET_1S.model", 3889380
    blob_sha1 = "2c2fff268e1dc50334c71ffbdd7792f6d21ed816"


ALL_ENGINES = (Egret1Engine, Egret1eEngine, Egret1tEngine, Egret1mEngine, Egret1sEngine)
