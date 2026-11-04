import numpy as np
import pytest

from iqmolx_engines.base import Context, EngineError, resolve_options
from iqmolx_engines.engines.demo import DemoEngine, RADII

WATER_Z = [8, 1, 1]
WATER_XYZ = np.array([[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]])


def make(**opts):
    eng = DemoEngine()
    eng.load(resolve_options(DemoEngine, opts, "load"), Context())
    return eng


def run(eng, xyz, want, charge=0, **opts):
    return eng.compute(WATER_Z, xyz, charge, 1, want,
                       resolve_options(DemoEngine, opts, "compute"), Context())


def test_forces_match_finite_difference_of_energy():
    eng = make()
    f = run(eng, WATER_XYZ, ["energy", "forces"])["forces"]
    h = 1e-5
    for a in range(3):
        for k in range(3):
            p, m = WATER_XYZ.copy(), WATER_XYZ.copy()
            p[a, k] += h
            m[a, k] -= h
            fd = -(run(eng, p, ["energy"])["energy"] - run(eng, m, ["energy"])["energy"]) / (2 * h)
            assert f[a, k] == pytest.approx(fd, abs=1e-6)


def test_hessian_matches_finite_difference_of_forces():
    eng = make()
    hess = run(eng, WATER_XYZ, ["energy", "hessian"])["hessian"]
    assert hess.shape == (9, 9)
    assert np.allclose(hess, hess.T)
    h = 1e-5
    fd = np.zeros((9, 9))
    for a in range(3):
        for k in range(3):
            p, m = WATER_XYZ.copy(), WATER_XYZ.copy()
            p[a, k] += h
            m[a, k] -= h
            fd[3 * a + k] = -(run(eng, p, ["forces"])["forces"]
                              - run(eng, m, ["forces"])["forces"]).ravel() / (2 * h)
    assert np.allclose(hess, fd, atol=1e-5)


def test_translation_invariance_and_zero_net_force():
    eng = make()
    r1 = run(eng, WATER_XYZ, ["energy", "forces"])
    r2 = run(eng, WATER_XYZ + [1.0, -2.0, 3.0], ["energy", "forces"])
    assert r1["energy"] == pytest.approx(r2["energy"])
    assert np.allclose(r1["forces"].sum(axis=0), 0.0, atol=1e-12)


def test_diatomic_minimum_at_sum_of_radii():
    eng = DemoEngine()
    eng.load(resolve_options(DemoEngine, {}, "load"), Context())
    r0 = RADII[1] * 2
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, r0]])
    res = eng.compute([1, 1], xyz, 0, 1, ["energy", "forces"],
                      resolve_options(DemoEngine, {}, "compute"), Context())
    assert res["energy"] == pytest.approx(-4.0)
    assert np.allclose(res["forces"], 0.0, atol=1e-12)


def test_dipole_and_charges():
    eng = make()
    res = run(eng, WATER_XYZ, ["dipole", "charges"], charge=1)
    assert res["charges"].sum() == pytest.approx(1.0)
    assert res["dipole"] == pytest.approx(res["charges"] @ WATER_XYZ)


def test_validation():
    with pytest.raises(EngineError) as e:
        DemoEngine.validate_system([92], 0, 1, ["energy"])
    assert e.value.code == "unsupported_element"
    with pytest.raises(EngineError) as e:
        DemoEngine.validate_system([1], 0, 1, ["spin_density"])
    assert e.value.code == "unsupported_property"


def test_option_checks():
    with pytest.raises(EngineError) as e:
        resolve_options(DemoEngine, {"well_depth": -1.0}, "load")
    assert e.value.code == "bad_option"
    with pytest.raises(EngineError) as e:
        resolve_options(DemoEngine, {"nonsense": 1}, "load")
    assert e.value.code == "bad_option"
    assert resolve_options(DemoEngine, {}, "load") == {"well_depth": 4.0, "alpha": 1.8}


def test_required_option_rules():
    from iqmolx_engines.base import Engine, option

    class Needy(Engine):
        name = "needy"
        options = (option("path", "Model path", "string", "", required=True),
                   option("level", "Level", "int", 0, required=True))

    with pytest.raises(EngineError) as e:
        resolve_options(Needy, {}, "load")
    assert e.value.code == "bad_option" and "Model path" in e.value.message
    assert resolve_options(Needy, {"path": "x"}, "load") == {"path": "x", "level": 0}


def test_options_that_are_shown_only_for_some_values():
    from iqmolx_engines.base import Engine, option, option_is_shown

    class Picky(Engine):
        name = "picky"
        options = (option("method", "Method", "choice", "a", choices=["a", "b", "other"]),
                   option("custom", "Custom name", "string", "", required=True, show_if={"method": "other"}),
                   option("grid", "Grid", "int", 1, show_if={"method": ["b", "other"]}))

    opts = {o["name"]: o for o in Picky.options}
    assert option_is_shown(opts["custom"], {"method": "other"}) and not option_is_shown(opts["custom"], {"method": "a"})
    assert option_is_shown(opts["grid"], {"method": "b"}) and option_is_shown(opts["grid"], {"method": "other"})
    assert not option_is_shown(opts["grid"], {"method": "a"})
    assert resolve_options(Picky, {}, "load")["custom"] == ""                  # hidden: not required
    assert resolve_options(Picky, {"method": "b"}, "load")["method"] == "b"
    for given in ({"method": "other"}, {"method": "other", "custom": "  "}):  # shown: required
        with pytest.raises(EngineError) as e:
            resolve_options(Picky, given, "load")
        assert e.value.code == "bad_option" and "Custom name" in e.value.message
    assert resolve_options(Picky, {"method": "other", "custom": "x"}, "load")["custom"] == "x"
