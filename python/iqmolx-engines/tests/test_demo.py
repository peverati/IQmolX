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
