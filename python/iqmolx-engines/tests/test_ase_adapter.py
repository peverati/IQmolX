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
