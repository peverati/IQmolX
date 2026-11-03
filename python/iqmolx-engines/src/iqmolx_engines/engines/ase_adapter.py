# Copyright (C) 2026 Roberto Peverati
# Developed with the assistance of Claude (Anthropic).
#
# This file is part of IQMolX, a fork of IQmol (http://iqmol.org). It is free
# software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version. It is
# distributed WITHOUT ANY WARRANTY; see <http://www.gnu.org/licenses/>.

"""Generic adapter: any ASE calculator as an engine. ASE is imported only when needed, so
the helper runs without it.

Charge and spin reach the calculator in two ways, because calculators differ:
* through `atoms.info` ("charge" and "spin" = multiplicity), the convention used by fairchem and
  AIMNet2;
* as the calculator's own parameters, for the ones it has: `charge`, and `multiplicity` or `mult`
  (TBLite, MOPAC, Q-Chem declare them; for ORCA or Gaussian they are there when the arguments
  contain them).  The values of the window replace those of the arguments at every calculation.
A calculator that reads neither would silently use a neutral singlet."""

import importlib
import json

import numpy as np

from ..base import Engine, EngineError, option


def _ase_missing():
    try:
        import ase  # noqa: F401
    except ImportError:
        return "ASE is not installed (pip install ase)"
    return ""


def _find_class(spec):
    """'package.module:Class' (or 'package.module.Class') to the class object."""
    module, sep, name = spec.partition(":")
    if not sep:
        module, _, name = spec.rpartition(".")
    if not module or not name:
        raise EngineError("bad_option", "Calculator must look like 'package.module:Class'")
    try:
        return getattr(importlib.import_module(module), name)
    except (ImportError, AttributeError) as exc:
        raise EngineError("engine_error", "Cannot import %s: %s" % (spec, exc)) from None


class AseEngine(Engine):
    name = "ase"
    title = "ASE calculator"
    description = ("Wraps any calculator of the Atomic Simulation Environment, for example "
                   "'ase.calculators.emt:EMT'. Constructor arguments are given as JSON. The charge "
                   "and spin of the window are passed on to calculators that have such arguments.")
    license = "depends on the calculator"
    elements = None
    charge = True
    multiplicity = True
    properties = ("energy", "forces")
    requirements = ("ase",)
    options = (
        option("calculator", "Calculator class", "string", "", required=True,
               placeholder="ase.calculators.emt:EMT",
               description="Python path 'package.module:Class'."),
        option("kwargs", "Constructor arguments (JSON)", "string", "{}",
               description='For example {"model_name": "x"}.'),
    )

    @classmethod
    def unavailable_reason(cls):
        return _ase_missing()

    def load(self, options, ctx):
        if not options["calculator"]:
            raise EngineError("bad_option", "Give the calculator class to use")
        try:
            kwargs = json.loads(options["kwargs"] or "{}")
        except ValueError:
            raise EngineError("bad_option", "Constructor arguments are not valid JSON") from None
        if not isinstance(kwargs, dict):
            raise EngineError("bad_option", "Constructor arguments must be a JSON object")
        cls = _find_class(options["calculator"])
        self.calc = cls(**kwargs)
        return {"calculator": options["calculator"]}

    def compute(self, numbers, coordinates, charge, multiplicity, want, options, ctx):
        from ase import Atoms
        atoms = Atoms(numbers=numbers, positions=coordinates)
        atoms.info["charge"] = charge
        atoms.info["spin"] = multiplicity
        self.pass_charge_and_spin(charge, multiplicity)
        atoms.calc = self.calc
        result = {"energy": atoms.get_potential_energy()}
        if "forces" in want:
            result["forces"] = np.asarray(atoms.get_forces())
        return result

    def pass_charge_and_spin(self, charge, multiplicity):
        """Set the calculator's charge and multiplicity parameters, if it has any."""
        names = set(getattr(self.calc, "parameters", None) or {})
        names |= set(getattr(type(self.calc), "default_parameters", None) or {})
        wanted = {}
        if "charge" in names:
            wanted["charge"] = charge
        for name in ("multiplicity", "mult"):
            if name in names:
                wanted[name] = multiplicity
        current = getattr(self.calc, "parameters", None) or {}
        changed = {k: v for k, v in wanted.items() if current.get(k) != v}
        if changed:
            self.calc.set(**changed)
            # ASE's base set() keeps old results; a calculator that does not reset itself would
            # otherwise answer with the numbers of the previous charge
            if hasattr(self.calc, "reset"):
                self.calc.reset()
