# Copyright (C) 2026 Roberto Peverati
# Developed with the assistance of Claude (Anthropic).
#
# This file is part of IQMolX, a fork of IQmol (http://iqmol.org). It is free
# software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version. It is
# distributed WITHOUT ANY WARRANTY; see <http://www.gnu.org/licenses/>.

"""TorchANI engine: the ANI-2x neural network potential (Smith, Isayev, Roitberg et al.),
through the TorchANI package (MIT license).  It covers H, C, N, O, F, S and Cl in neutral,
closed-shell molecules.  PyTorch and TorchANI are imported only when the model is loaded, so
the helper starts without them.

The parameters are not part of the TorchANI package.  TorchANI downloads ani2x_state_dict.pt
the first time the model is used, from the public Hugging Face repository
roitberg-group/ani2x (MIT license, no login needed), and keeps it in its own cache."""

import importlib.util

import numpy as np

from ..base import Engine, EngineError, option

HARTREE_TO_EV = 27.211386245988

# The elements of ANI-2x: H, C, N, O, F, S, Cl
ELEMENTS = (1, 6, 7, 8, 9, 16, 17)


class TorchAniEngine(Engine):
    name = "torchani"
    title = "ANI-2x (TorchANI)"
    description = ("ANI-2x neural network potential, trained on wB97X/6-31G(d) energies and "
                   "forces. Neutral closed-shell molecules of H, C, N, O, F, S and Cl only.")
    license = "MIT (TorchANI and the ANI-2x parameters)"
    citation = ("Devereux et al., J. Chem. Theory Comput. 2020, 16, 4192; "
                "Gao et al., J. Chem. Inf. Model. 2020, 60, 3408")
    elements = ELEMENTS
    charge = False
    multiplicity = False
    properties = ("energy", "forces", "hessian")
    requires_download = True
    requirements = ("torch", "torchani")
    options = (
        option("model_index", "Network (-1 for the ensemble)", "int", -1, min=-1, max=7,
               description="ANI-2x is an ensemble of 8 networks. -1 averages all of them "
                           "(the standard, slower); 0 to 7 use a single network."),
        option("precision", "Precision", "choice", "float64", choices=["float64", "float32"],
               description="float64 gives smoother forces for optimizations and frequencies."),
    )

    @classmethod
    def unavailable_reason(cls):
        missing = [m for m in ("torch", "torchani") if importlib.util.find_spec(m) is None]
        if missing:
            return "%s not installed (pip install torch torchani)" % " and ".join(missing)
        return ""

    def load(self, options, ctx):
        import torch
        import torchani

        ctx.progress(None, "Loading ANI-2x (the parameters are downloaded the first time)")
        ctx.check_cancel()
        self.torch = torch
        self.dtype = torch.float64 if options["precision"] == "float64" else torch.float32
        index = options["model_index"]
        try:
            model = torchani.models.ANI2x(periodic_table_index=True, device="cpu",
                                          dtype=self.dtype)
        except Exception as exc:
            raise EngineError("engine_error", "Cannot load the ANI-2x model: %s" % exc,
                              detail="The parameters come from huggingface.co/roitberg-group/"
                                     "ani2x; check the internet connection.") from None
        if index >= 0:
            model = model[index]
        self.model = model
        return {"device": "cpu", "model": "ANI-2x" + (" ensemble" if index < 0 else " #%d" % index),
                "torch": torch.__version__, "torchani": getattr(torchani, "__version__", "")}

    def energy_of(self, species, coords):
        """Energy in Hartree of one molecule (coords: 1 x N x 3 tensor)."""
        return self.model((species, coords)).energies.sum()

    def compute(self, numbers, coordinates, charge, multiplicity, want, options, ctx):
        torch = self.torch
        ctx.check_cancel()
        n = len(numbers)
        species = torch.tensor(np.asarray(numbers, dtype=np.int64)).unsqueeze(0)
        xyz = torch.tensor(np.asarray(coordinates), dtype=self.dtype).unsqueeze(0)
        xyz.requires_grad_(True)

        energy = self.energy_of(species, xyz)
        result = {"energy": energy.item() * HARTREE_TO_EV}
        if "forces" in want:
            grad, = torch.autograd.grad(energy, xyz)
            result["forces"] = -grad.detach().numpy().reshape(n, 3) * HARTREE_TO_EV
        if "hessian" in want:
            ctx.check_cancel()
            flat = torch.tensor(np.asarray(coordinates), dtype=self.dtype).reshape(-1)
            hess = torch.autograd.functional.hessian(
                lambda x: self.energy_of(species, x.reshape(1, n, 3)), flat)
            result["hessian"] = hess.detach().numpy() * HARTREE_TO_EV
        return result
