import os

import numpy as np
import pytest

from iqmolx_engines.base import Context, EngineError, resolve_options
from iqmolx_engines.engines.torchani_engine import ALL_ENGINES, TorchAniEngine

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
    h = 1e-5     # the finite-difference error only vanishes for small steps (float64)
    fd = np.zeros((9, 9))
    for a in range(3):
        for k in range(3):
            p, m = WATER_XYZ.copy(), WATER_XYZ.copy()
            p[a, k] += h
            m[a, k] -= h
            fd[3 * a + k] = -(run(engine, p, ["forces"])["forces"]
                              - run(engine, m, ["forces"])["forces"]).ravel() / (2 * h)
    assert np.allclose(hess, fd, atol=1e-4)


@needs_model
def test_water_is_a_sensible_molecule(engine):
    # not a reference value, only a check that the energy is that of water and not noise:
    # translation invariance and nearly vanishing net force
    r1 = run(engine, WATER_XYZ, ["energy", "forces"])
    r2 = run(engine, WATER_XYZ + [3.0, -1.0, 2.0], ["energy", "forces"])
    assert r1["energy"] == pytest.approx(r2["energy"], abs=1e-6)
    assert np.allclose(r1["forces"].sum(axis=0), 0.0, atol=1e-3)


# --- all the ANI models ------------------------------------------------------------------------

def test_every_model_describes_itself_properly():
    names = [c.name for c in ALL_ENGINES]
    assert len(set(names)) == len(names) and "torchani" in names       # unique; ANI-2x keeps its name
    for cls in ALL_ENGINES:
        d = cls.describe()
        assert d["name"] and d["title"] and d["description"] and d["citation"]
        assert d["elements"] in ([1, 6, 7, 8], [1, 6, 7, 8, 9, 16, 17])
        assert d["charge"] is False and d["multiplicity"] is False
        assert d["requires_download"] is True and d["requirements"] == ["torch", "torchani"]
        assert {"energy", "forces", "hessian"} <= set(d["properties"])
        assert "MIT" in d["license"]
        assert [o["name"] for o in d["options"]][:2] == ["model_index", "precision"]


def test_elements_follow_the_models():
    by = {c.name: c for c in ALL_ENGINES}
    for name in ("ani1x", "ani1ccx", "ani1xnr"):
        assert by[name].elements == (1, 6, 7, 8)                       # HCNO only
        with pytest.raises(EngineError) as e:
            by[name].validate_system([9, 1], 0, 1, ["energy"])         # fluorine
        assert e.value.code == "unsupported_element"
    for name in ("torchani", "ani2xr", "ani2dr", "anir2s", "animbis"):
        assert by[name].elements == (1, 6, 7, 8, 9, 16, 17)


def test_the_special_models():
    by = {c.name: c for c in ALL_ENGINES}
    assert {"charges", "dipole"} <= set(by["animbis"].properties)
    solvent = [o for o in by["anir2s"].options if o["name"] == "solvent"][0]
    assert solvent["choices"] == ["vacuum", "water", "chcl3", "ch3cn"] and solvent["default"] == "vacuum"
    assert by["anir2s"]().model_arguments({"solvent": "vacuum"}) == {"solvent": None}
    assert by["anir2s"]().model_arguments({"solvent": "water"}) == {"solvent": "water"}


# The models themselves need torch, torchani and their parameters (10 to 55 MB each, downloaded the
# first time).  IQMOLX_TEST_MODELS names the engines to run, e.g. IQMOLX_TEST_MODELS=torchani,ani1x
# (or "all"); it needs IQMOLX_TEST_DOWNLOADS=1 as well.
def _models_asked_for():
    value = os.environ.get("IQMOLX_TEST_MODELS", "")
    if os.environ.get("IQMOLX_TEST_DOWNLOADS") != "1" or TorchAniEngine.unavailable_reason():
        return []
    names = [c.name for c in ALL_ENGINES] if value == "all" else [v for v in value.split(",") if v]
    return [c for c in ALL_ENGINES if c.name in names]


@pytest.mark.parametrize("cls", _models_asked_for(), ids=lambda c: c.name)
def test_each_model_gives_sensible_water(cls):
    eng = cls()
    eng.load(resolve_options(cls, {}, "load"), Context())
    want = ["energy", "forces", "hessian"] + [p for p in ("charges", "dipole") if p in cls.properties]
    res = eng.compute(WATER_Z, WATER_XYZ, 0, 1, want, {}, Context())
    # forces are minus the gradient of the energy
    h = 1e-4
    for a, k in ((0, 2), (1, 1), (2, 1)):
        p, m = WATER_XYZ.copy(), WATER_XYZ.copy()
        p[a, k] += h
        m[a, k] -= h
        ep = eng.compute(WATER_Z, p, 0, 1, ["energy"], {}, Context())["energy"]
        em = eng.compute(WATER_Z, m, 0, 1, ["energy"], {}, Context())["energy"]
        assert res["forces"][a, k] == pytest.approx(-(ep - em) / (2 * h), abs=2e-4)
    assert res["hessian"].shape == (9, 9) and np.allclose(res["hessian"], res["hessian"].T, atol=1e-6)
    assert np.allclose(res["forces"].sum(axis=0), 0.0, atol=2e-3)
    assert -2300 < res["energy"] < -2000                 # water, in eV (about -76 Hartree)
    if "charges" in res:
        assert res["charges"].sum() == pytest.approx(0.0, abs=0.15) and res["charges"][0] < 0
        assert res["dipole"].shape == (3,)


def test_the_solvent_model_changes_the_energy():
    if "anir2s" not in [c.name for c in _models_asked_for()]:
        pytest.skip("set IQMOLX_TEST_DOWNLOADS=1 and IQMOLX_TEST_MODELS=anir2s")
    from iqmolx_engines.engines.torchani_engine import AniR2sEngine
    energies = {}
    for solvent in ("vacuum", "water"):
        eng = AniR2sEngine()
        eng.load(resolve_options(AniR2sEngine, {"solvent": solvent}, "load"), Context())
        energies[solvent] = eng.compute(WATER_Z, WATER_XYZ, 0, 1, ["energy"], {}, Context())["energy"]
    assert energies["water"] < energies["vacuum"]        # solvation stabilizes water
