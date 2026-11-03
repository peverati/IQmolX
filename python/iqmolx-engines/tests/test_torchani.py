import os

import numpy as np
import pytest

from iqmolx_engines.base import Context, EngineError, resolve_options
from iqmolx_engines.engines.torchani_engine import TorchAniEngine

WATER_Z = [8, 1, 1]
WATER_XYZ = np.array([[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]])


def test_description_and_capabilities():
    d = TorchAniEngine.describe()
    assert d["name"] == "torchani" and d["requires_download"] is True
    assert d["elements"] == [1, 6, 7, 8, 9, 16, 17]
    assert d["charge"] is False and d["multiplicity"] is False
    assert {"energy", "forces", "hessian"} <= set(d["properties"])
    assert d["available"] == (d["unavailable_reason"] == "")
    assert resolve_options(TorchAniEngine, {}, "load") == {"model_index": -1,
                                                           "precision": "float64"}


def test_unavailable_reason_names_what_is_missing():
    reason = TorchAniEngine.unavailable_reason()
    try:
        import torch, torchani  # noqa: F401
        assert reason == ""
    except ImportError:
        assert "pip install" in reason


def test_what_ani2x_cannot_do_is_refused_before_the_model_runs():
    with pytest.raises(EngineError) as e:
        TorchAniEngine.validate_system([26, 1], 0, 1, ["energy"])       # iron
    assert e.value.code == "unsupported_element"
    with pytest.raises(EngineError) as e:
        TorchAniEngine.validate_system(WATER_Z, 1, 1, ["energy"])
    assert e.value.code == "unsupported_charge"
    with pytest.raises(EngineError) as e:
        TorchAniEngine.validate_system(WATER_Z, 0, 3, ["energy"])
    assert e.value.code == "unsupported_multiplicity"


# The model itself needs torch, torchani and a download of the parameters (a few MB from
# Hugging Face), so it only runs when asked for: IQMOLX_TEST_DOWNLOADS=1 pytest
needs_model = pytest.mark.skipif(
    os.environ.get("IQMOLX_TEST_DOWNLOADS") != "1" or TorchAniEngine.unavailable_reason() != "",
    reason="set IQMOLX_TEST_DOWNLOADS=1 with torch and torchani installed")


@pytest.fixture(scope="module")
def engine():
    eng = TorchAniEngine()
    eng.load(resolve_options(TorchAniEngine, {}, "load"), Context())
    return eng


def run(eng, xyz, want):
    return eng.compute(WATER_Z, xyz, 0, 1, want, {}, Context())


@needs_model
def test_forces_match_finite_differences(engine):
    f = run(engine, WATER_XYZ, ["energy", "forces"])["forces"]
    h = 1e-4
    for a in range(3):
        for k in range(3):
            p, m = WATER_XYZ.copy(), WATER_XYZ.copy()
            p[a, k] += h
            m[a, k] -= h
            fd = -(run(engine, p, ["energy"])["energy"] - run(engine, m, ["energy"])["energy"]) / (2 * h)
            assert f[a, k] == pytest.approx(fd, abs=1e-4)


@needs_model
def test_hessian_matches_finite_differences_of_forces(engine):
    hess = run(engine, WATER_XYZ, ["energy", "hessian"])["hessian"]
    assert hess.shape == (9, 9) and np.allclose(hess, hess.T, atol=1e-6)
    h = 1e-4
    fd = np.zeros((9, 9))
    for a in range(3):
        for k in range(3):
            p, m = WATER_XYZ.copy(), WATER_XYZ.copy()
            p[a, k] += h
            m[a, k] -= h
            fd[3 * a + k] = -(run(engine, p, ["forces"])["forces"]
                              - run(engine, m, ["forces"])["forces"]).ravel() / (2 * h)
    assert np.allclose(hess, fd, atol=1e-2)


@needs_model
def test_water_is_a_sensible_molecule(engine):
    # not a reference value, only a check that the energy is that of water and not noise:
    # translation invariance and nearly vanishing net force
    r1 = run(engine, WATER_XYZ, ["energy", "forces"])
    r2 = run(engine, WATER_XYZ + [3.0, -1.0, 2.0], ["energy", "forces"])
    assert r1["energy"] == pytest.approx(r2["energy"], abs=1e-6)
    assert np.allclose(r1["forces"].sum(axis=0), 0.0, atol=1e-3)
