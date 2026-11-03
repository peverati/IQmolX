import numpy as np
import pytest

from iqmolx_engines.base import Context, EngineError, resolve_options
from iqmolx_engines.engines.ase_adapter import AseEngine, _find_class


def test_description_always_works():
    d = AseEngine.describe()
    assert d["name"] == "ase" and d["elements"] is None
    assert d["available"] == (AseEngine.unavailable_reason() == "")


def test_class_path_forms():
    assert _find_class("collections:OrderedDict").__name__ == "OrderedDict"
    assert _find_class("collections.OrderedDict").__name__ == "OrderedDict"
    with pytest.raises(EngineError) as e:
        _find_class("no_such_module_xyz:Thing")
    assert e.value.code == "engine_error"
    with pytest.raises(EngineError) as e:
        _find_class("justaname")
    assert e.value.code == "bad_option"


def test_load_validates_options():
    eng = AseEngine()
    with pytest.raises(EngineError) as e:
        eng.load(resolve_options(AseEngine, {}, "load"), Context())
    assert e.value.code == "bad_option"
    with pytest.raises(EngineError) as e:
        eng.load(resolve_options(AseEngine, {"calculator": "collections:OrderedDict",
                                             "kwargs": "[1]"}, "load"), Context())
    assert e.value.code == "bad_option"


def test_emt_through_adapter():
    pytest.importorskip("ase")
    eng = AseEngine()
    eng.load(resolve_options(AseEngine, {"calculator": "ase.calculators.emt:EMT"}, "load"),
             Context())
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 2.5]])
    res = eng.compute([29, 29], xyz, 0, 1, ["energy", "forces"], {}, Context())
    assert np.isfinite(res["energy"])
    assert res["forces"].shape == (2, 3)
    assert np.allclose(res["forces"].sum(axis=0), 0.0, atol=1e-10)


# --- charge and spin reach calculators that take them as arguments ---------------------------

def _probe_class():
    pytest.importorskip("ase")
    from ase.calculators.calculator import Calculator

    class ProbeCalculator(Calculator):
        """Energy = 100 * charge + multiplicity, so the test sees what the calculator received."""
        implemented_properties = ["energy", "forces"]
        default_parameters = {"charge": None, "multiplicity": None}

        def calculate(self, atoms=None, properties=("energy",), system_changes=()):
            Calculator.calculate(self, atoms, properties, system_changes)
            q, m = self.parameters.charge, self.parameters.multiplicity
            self.results = {"energy": 100.0 * (q or 0) + (m or 0), "forces": np.zeros((len(atoms), 3))}
    return ProbeCalculator


def test_charge_and_multiplicity_are_set_on_the_calculator(monkeypatch):
    probe = _probe_class()
    import sys, types
    module = types.ModuleType("tests_probe")
    module.ProbeCalculator = probe
    monkeypatch.setitem(sys.modules, "tests_probe", module)
    eng = AseEngine()
    eng.load(resolve_options(AseEngine, {"calculator": "tests_probe:ProbeCalculator"}, "load"), Context())
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    assert eng.compute([1, 1], xyz, 0, 1, ["energy"], {}, Context())["energy"] == pytest.approx(1.0)
    assert eng.compute([1, 1], xyz, 1, 2, ["energy"], {}, Context())["energy"] == pytest.approx(102.0)
    assert eng.compute([1, 1], xyz, -1, 3, ["energy"], {}, Context())["energy"] == pytest.approx(-97.0)


def test_tblite_gets_the_charge():
    pytest.importorskip("ase")
    pytest.importorskip("tblite")
    eng = AseEngine()
    eng.load(resolve_options(AseEngine, {"calculator": "tblite.ase:TBLite",
                                         "kwargs": '{"method": "GFN2-xTB", "verbosity": 0}'},
                             "load"), Context())
    water = [[0, 0, 0.119], [0, 0.763, -0.477], [0, -0.763, -0.477]]
    neutral = eng.compute([8, 1, 1], np.array(water), 0, 1, ["energy"], {}, Context())["energy"]
    cation = eng.compute([8, 1, 1], np.array(water), 1, 2, ["energy"], {}, Context())["energy"]
    assert neutral == pytest.approx(-137.968, abs=2e-3)       # GFN2-xTB, the value in the tblite docs
    assert cation - neutral > 5.0                              # it was silently neutral before
