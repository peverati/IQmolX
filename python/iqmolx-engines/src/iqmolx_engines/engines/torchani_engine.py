# Copyright (C) 2026 Roberto Peverati
# Developed with the assistance of Claude (Anthropic).
#
# This file is part of IQMolX, a fork of IQmol (http://iqmol.org). It is free
# software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version. It is
# distributed WITHOUT ANY WARRANTY; see <http://www.gnu.org/licenses/>.
"""TorchANI engines: the ANI family of neural network potentials (Smith, Isayev, Roitberg and
others) through the TorchANI package (MIT license).  One engine per model, so that each one says
exactly which elements, level of theory and license it has.  PyTorch and TorchANI are imported
only when a model is loaded, so the helper starts without them.

The parameters are not part of the TorchANI package.  TorchANI downloads them the first time a
model is used, from public repositories of the Roitberg group on Hugging Face (roitberg-group/*,
MIT license on their model cards, no login needed), and keeps them in its own cache."""

import importlib.util

import numpy as np

from ..base import Engine, EngineError, option

HARTREE_TO_EV = 27.211386245988

HCNO = (1, 6, 7, 8)                      # H, C, N, O
ELEMENTS_2X = (1, 6, 7, 8, 9, 16, 17)    # H, C, N, O, F, S, Cl

GAO = "Gao et al., J. Chem. Inf. Model. 60, 3408 (2020) (TorchANI)"


class TorchAniBase(Engine):
    """What all the ANI models share.  A model sets `name`, `title`, `description`, `elements`,
    `citation` and `model_function`, the function of torchani.models that builds it."""

    model_function = ""
    license = "MIT (TorchANI and the parameters)"
    charge = False
    multiplicity = False
    properties = ("energy", "forces", "hessian")
    requires_download = True
    requirements = ("torch", "torchani")
    options = (
        option("model_index", "Network (-1 for the ensemble)", "int", -1, min=-1, max=7,
               description="The models are ensembles of 8 networks. -1 averages all of them "
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

    def model_arguments(self, options):
        """Extra arguments of the torchani.models function (a model may have settings)."""
        return {}

    def load(self, options, ctx):
        import torch
        import torchani

        ctx.progress(None, "Loading %s (the parameters are downloaded the first time)" % self.title)
        ctx.check_cancel()
        self.torch = torch
        self.dtype = torch.float64 if options["precision"] == "float64" else torch.float32
        index = options["model_index"]
        build = getattr(torchani.models, self.model_function)
        try:
            model = build(periodic_table_index=True, device="cpu", dtype=self.dtype,
                          **self.model_arguments(options))
        except Exception as exc:
            raise EngineError("engine_error", "Cannot load the %s model: %s" % (self.title, exc),
                              detail="The parameters come from huggingface.co/roitberg-group; "
                                     "check the internet connection.") from None
        if index >= 0:
            model = model[index]
        self.model = model
        return {"device": "cpu", "model": self.title + (" ensemble" if index < 0 else " #%d" % index),
                "torch": torch.__version__, "torchani": getattr(torchani, "__version__", "")}

    def evaluate(self, species, coords):
        """Energy in Hartree of one molecule (coords: 1 x N x 3 tensor) and, for a model that
        has them, the atomic charges (else None)."""
        out = self.model((species, coords))
        return out.energies.sum(), getattr(out, "atomic_charges", None)

    def compute(self, numbers, coordinates, charge, multiplicity, want, options, ctx):
        torch = self.torch
        ctx.check_cancel()
        n = len(numbers)
        species = torch.tensor(np.asarray(numbers, dtype=np.int64)).unsqueeze(0)
        xyz = torch.tensor(np.asarray(coordinates), dtype=self.dtype).unsqueeze(0)
        xyz.requires_grad_(True)

        energy, charges = self.evaluate(species, xyz)
        result = {"energy": energy.item() * HARTREE_TO_EV}
        if "forces" in want:
            grad, = torch.autograd.grad(energy, xyz, retain_graph=charges is not None)
            result["forces"] = -grad.detach().numpy().reshape(n, 3) * HARTREE_TO_EV
        if charges is not None:
            q = charges.detach().numpy().reshape(n)
            if "charges" in want:
                result["charges"] = q
            if "dipole" in want:
                result["dipole"] = q @ np.asarray(coordinates, dtype=float)    # e*Angstrom
        if "hessian" in want:
            ctx.check_cancel()
            flat = torch.tensor(np.asarray(coordinates), dtype=self.dtype).reshape(-1)
            hess = torch.autograd.functional.hessian(
                lambda x: self.evaluate(species, x.reshape(1, n, 3))[0], flat)
            result["hessian"] = hess.detach().numpy() * HARTREE_TO_EV
        return result


class TorchAniEngine(TorchAniBase):
    """ANI-2x.  Its name stays `torchani`, as saved preferences use it."""
    name = "torchani"
    title = "ANI-2x (TorchANI)"
    description = ("ANI-2x neural network potential, trained on wB97X/6-31G(d) energies and "
                   "forces. Neutral closed-shell molecules of H, C, N, O, F, S and Cl only. "
                   "Downloads about 55 MB the first time.")
    citation = ("Devereux et al., J. Chem. Theory Comput. 2020, 16, 4192; " + GAO)
    elements = ELEMENTS_2X
    model_function = "ANI2x"


class Ani1xEngine(TorchAniBase):
    name = "ani1x"
    title = "ANI-1x (TorchANI)"
    description = ("ANI-1x neural network potential, trained by active learning to "
                   "wB97X/6-31G(d). Neutral closed-shell molecules of H, C, N and O only. "
                   "Downloads about 10 MB the first time.")
    citation = "Smith et al., J. Chem. Phys. 148, 241733 (2018); " + GAO
    elements = HCNO
    model_function = "ANI1x"


class Ani1ccxEngine(TorchAniBase):
    name = "ani1ccx"
    title = "ANI-1ccx (TorchANI)"
    description = ("ANI-1ccx neural network potential, transfer-learned to CCSD(T)*/CBS accuracy "
                   "from ANI-1x. Neutral closed-shell molecules of H, C, N and O only. "
                   "Downloads about 10 MB the first time.")
    citation = "Smith et al., Nat. Commun. 10, 2903 (2019); " + GAO
    elements = HCNO
    model_function = "ANI1ccx"


class Ani1xnrEngine(TorchAniBase):
    name = "ani1xnr"
    title = "ANI-1xnr (TorchANI)"
    description = ("ANI-1xnr, a reactive neural network potential trained to BLYP/TZV2P. "
                   "Neutral closed-shell molecules of H, C, N and O only. Its energies are on a "
                   "different scale from the other ANI models (water is about -468 eV, not -2079 eV): "
                   "only energy differences mean something. Downloads about 21 MB the first time.")
    citation = "Zhang et al., Nat. Chem. (2024), doi:10.1038/s41557-023-01427-3; " + GAO
    elements = HCNO
    model_function = "ANI1xnr"


class Ani2xrEngine(TorchAniBase):
    name = "ani2xr"
    title = "ANI-2xr (TorchANI)"
    description = ("ANI-2xr: the ANI-2x data and wB97X/6-31G(d) level with an added repulsion "
                   "potential and a smoother potential energy surface. Neutral closed-shell "
                   "molecules of H, C, N, O, F, S and Cl only. Downloads about 55 MB the first time.")
    citation = "See the TorchANI repository for this model; " + GAO
    elements = ELEMENTS_2X
    model_function = "ANI2xr"


class Ani2drEngine(TorchAniBase):
    name = "ani2dr"
    title = "ANI-2dr (TorchANI)"
    description = ("ANI-2dr: trained to the B97-3c level on the ANI-2x data, with added "
                   "dispersion and repulsion potentials and a smoother potential energy surface. "
                   "Neutral closed-shell molecules of H, C, N, O, F, S and Cl only. "
                   "Downloads about 55 MB the first time.")
    citation = "See the TorchANI repository for this model; " + GAO
    elements = ELEMENTS_2X
    model_function = "ANI2dr"


class AniR2sEngine(TorchAniBase):
    name = "anir2s"
    title = "ANI-r2s (TorchANI, experimental)"
    description = ("ANI-r2s, an experimental model trained to r2SCAN-3c, in vacuum or with an "
                   "implicit solvent (SMD): water, chloroform or acetonitrile. Neutral closed-shell "
                   "molecules of H, C, N, O, F, S and Cl only. Downloads about 55 MB the first time "
                   "for each solvent.")
    citation = "See the TorchANI repository for this model; " + GAO
    elements = ELEMENTS_2X
    model_function = "ANIr2s"
    options = TorchAniBase.options + (
        option("solvent", "Solvent", "choice", "vacuum", choices=["vacuum", "water", "chcl3", "ch3cn"],
               description="vacuum, or the implicit solvent (SMD) the model was trained with: "
                           "water, chloroform (chcl3) or acetonitrile (ch3cn)."),
    )

    def model_arguments(self, options):
        return {"solvent": None if options["solvent"] == "vacuum" else options["solvent"]}


class AniMbisEngine(TorchAniBase):
    name = "animbis"
    title = "ANI-2x with MBIS charges (TorchANI, experimental)"
    description = ("An experimental ANI-2x that also predicts MBIS atomic charges. Neutral closed-shell "
                   "molecules of H, C, N, O, F, S and Cl only. The dipole is calculated from the "
                   "point charges. Downloads about 62 MB the first time.")
    citation = "See the TorchANI repository for this model; " + GAO
    elements = ELEMENTS_2X
    model_function = "ANImbis"
    properties = ("energy", "forces", "hessian", "charges", "dipole")


ALL_ENGINES = (TorchAniEngine, Ani1xEngine, Ani1ccxEngine, Ani1xnrEngine, Ani2xrEngine,
               Ani2drEngine, AniR2sEngine, AniMbisEngine)
