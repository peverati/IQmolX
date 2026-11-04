import numpy as np
import pytest

from iqmolx_engines.base import Cancelled, Context, EngineError, resolve_options
from iqmolx_engines.engines.pyscf_engine import METHODS, PySCFEngine

HARTREE = 27.211386245988
WATER_Z = [8, 1, 1]
WATER_XYZ = np.array([[0.0, 0.0, 0.119], [0.0, 0.763, -0.477], [0.0, -0.763, -0.477]])


def test_descriptor():
    d = PySCFEngine.describe()
    assert d["elements"] is None and d["charge"] is True and d["multiplicity"] is True
    assert {"energy", "forces", "hessian", "dipole", "charges"} <= set(d["properties"])
    assert d["requirements"] == ["pyscf"] and "Apache" in d["license"]
    basis = [o for o in d["options"] if o["name"] == "basis"][0]
    assert basis["required"] is True and basis["default"] == "def2-svp"
    method = [o for o in d["options"] if o["name"] == "method"][0]
    assert method["choices"][0] == "hf" and "wb97m-v" in method["choices"] and method["choices"][-1] == "other"
    other = [o for o in d["options"] if o["name"] == "method_other"][0]
    assert other["required"] is True and other["show_if"] == {"method": "other"}
    nlc = [o for o in d["options"] if o["name"] == "nlc_grid_level"][0]
    assert nlc["default"] == 1 and nlc["show_if"] == {"method": ["wb97m-v", "other"]}
    assert d["available"] == (d["unavailable_reason"] == "")


def test_the_name_of_an_other_functional_is_required_only_when_other_is_chosen():
    assert resolve_options(PySCFEngine, {}, "load")["method"] == "b3lyp"                # hidden, not required
    for given in ({"method": "other"}, {"method": "other", "method_other": "  "}):
        with pytest.raises(EngineError) as e:
            resolve_options(PySCFEngine, given, "load")
        assert e.value.code == "bad_option" and "Functional" in e.value.message
    assert resolve_options(PySCFEngine, {"method": "other", "method_other": "cam-b3lyp"}, "load")


def test_an_empty_basis_is_refused_by_the_helper_rules():
    with pytest.raises(EngineError) as e:
        resolve_options(PySCFEngine, {"basis": "  "}, "load")
    assert e.value.code == "bad_option" and "Basis set" in e.value.message


pyscf = pytest.importorskip("pyscf")


def make(**opts):
    eng = PySCFEngine()
    eng.load(resolve_options(PySCFEngine, opts, "load"), Context())
    return eng


def run(eng, z, xyz, want=("energy",), charge=0, mult=1):
    return eng.compute(z, np.asarray(xyz, float), charge, mult, list(want), {}, Context())


def test_hydrogen_molecule_matches_the_textbook_value():
    e = run(make(method="hf", basis="sto-3g", density_fit=False), [1, 1], [[0, 0, 0], [0, 0, 0.74]])["energy"]
    assert e / HARTREE == pytest.approx(-1.116759, abs=2e-6)           # RHF/STO-3G at 0.74 Angstrom


@pytest.mark.parametrize("method", ["hf", "b3lyp", "tpss", "r2scan"])
def test_forces_and_hessian_are_derivatives_of_the_energy(method):
    eng = make(method=method, basis="6-31g*", density_fit=False)
    res = run(eng, WATER_Z, WATER_XYZ, ["energy", "forces", "hessian"])
    h = 2e-3
    for a, k in ((0, 2), (1, 1), (2, 1)):
        p, m = WATER_XYZ.copy(), WATER_XYZ.copy()
        p[a, k] += h
        m[a, k] -= h
        fd = -(run(eng, WATER_Z, p)["energy"] - run(eng, WATER_Z, m)["energy"]) / (2 * h)
        assert res["forces"][a, k] == pytest.approx(fd, abs=2e-3)
    assert np.allclose(res["forces"].sum(axis=0), 0.0, atol=1e-3)
    hess = res["hessian"]
    assert hess.shape == (9, 9) and np.allclose(hess, hess.T, atol=1e-4)
    for a, k in ((0, 2), (1, 1)):
        p, m = WATER_XYZ.copy(), WATER_XYZ.copy()
        p[a, k] += h
        m[a, k] -= h
        fd = -(run(eng, WATER_Z, p, ["forces"])["forces"] - run(eng, WATER_Z, m, ["forces"])["forces"]).ravel() / (2 * h)
        assert np.allclose(hess[3 * a + k], fd, atol=8e-2)        # was off by 0.6 to 6 for the meta-GGAs


def test_translation_invariance_dipole_and_charges():
    eng = make(method="hf", basis="def2-svp")
    r1 = run(eng, WATER_Z, WATER_XYZ, ["energy", "dipole", "charges"])
    r2 = run(eng, WATER_Z, WATER_XYZ + [1.5, -2.0, 0.7], ["energy", "dipole", "charges"])
    assert r1["energy"] == pytest.approx(r2["energy"], abs=1e-6)
    assert r1["charges"].sum() == pytest.approx(0.0, abs=1e-6) and r1["charges"][0] < 0
    assert r1["dipole"] == pytest.approx(r2["dipole"], abs=1e-4) and np.linalg.norm(r1["dipole"]) > 0.3


def test_the_ionization_energy_of_water_is_about_right():
    eng = make(method="b3lyp", basis="def2-svp")
    neutral = run(eng, WATER_Z, WATER_XYZ)["energy"]
    cation = run(eng, WATER_Z, WATER_XYZ, charge=1, mult=2)["energy"]      # open shell: unrestricted
    assert 11.5 < cation - neutral < 13.5                                     # eV; the experiment gives 12.6


def test_open_shell_molecule_and_the_reuse_of_the_density():
    eng = make(method="b3lyp", basis="def2-svp")
    oh = [[0, 0, 0], [0, 0, 0.97]]
    first = run(eng, [8, 1], oh, mult=2)["energy"]
    again = run(eng, [8, 1], oh, mult=2)["energy"]                          # starts from the last density
    fresh = run(make(method="b3lyp", basis="def2-svp"), [8, 1], oh, mult=2)["energy"]
    assert first == pytest.approx(again, abs=1e-6) == pytest.approx(fresh, abs=1e-6)
    assert first / HARTREE == pytest.approx(-75.667, abs=2e-3)


def test_density_fitting_is_close():
    exact = run(make(method="b3lyp", basis="def2-svp", density_fit=False), WATER_Z, WATER_XYZ)["energy"]
    fitted = run(make(method="b3lyp", basis="def2-svp", density_fit=True), WATER_Z, WATER_XYZ)["energy"]
    assert abs(exact - fitted) / HARTREE < 5e-4


def test_impossible_spin_unknown_basis_and_cancel_are_errors_not_crashes():
    eng = make(method="hf", basis="sto-3g")
    with pytest.raises(EngineError) as e:
        run(eng, WATER_Z, WATER_XYZ, mult=2)                                # 10 electrons, a doublet
    assert e.value.code == "unsupported_multiplicity" and "10 electrons" in e.value.message
    with pytest.raises(EngineError) as e:
        run(make(basis="no-such-basis"), WATER_Z, WATER_XYZ)
    assert e.value.code == "engine_error" and "no-such-basis" in e.value.message
    with pytest.raises(Cancelled):
        eng.compute(WATER_Z, WATER_XYZ, 0, 1, ["energy"], {}, Context(cancelled=lambda: True))


def test_a_heavy_atom_with_a_pseudopotential_works():
    eng = make(method="hf", basis="def2-svp")
    res = run(eng, [6, 53, 1, 1, 1], [[0, 0, 0], [0, 0, 2.14], [1.03, 0, -0.36], [-0.51, 0.89, -0.36],
                                       [-0.51, -0.89, -0.36]], ["energy", "forces"])
    assert np.isfinite(res["energy"]) and res["forces"].shape == (5, 3)


def test_the_numerical_hessian_agrees_with_the_analytic_one():
    eng = make(method="b3lyp", basis="6-31g*")
    analytic = run(eng, WATER_Z, WATER_XYZ, ["hessian"])["hessian"]
    numerical = eng.numerical_hessian(WATER_Z, WATER_XYZ, 0, 1, Context())
    assert np.abs(analytic - numerical).max() < 0.1
    assert not eng.meta_gga and make(method="tpss").meta_gga and make(method="wb97m-v").meta_gga


def test_wb97m_v_in_the_menu_gives_the_expected_energy_and_forces():
    eng = make(method="wb97m-v", basis="def2-svp")            # nonlocal grid level 1 (the default of this setting)
    res = run(eng, WATER_Z, WATER_XYZ, ["energy", "forces"])
    assert res["energy"] / HARTREE == pytest.approx(-76.32554, abs=2e-6)     # PySCF at its own finer VV10 grid
    h = 2e-3
    for a, k in ((0, 2), (1, 1)):
        p, m = WATER_XYZ.copy(), WATER_XYZ.copy()
        p[a, k] += h
        m[a, k] -= h
        fd = -(run(eng, WATER_Z, p)["energy"] - run(eng, WATER_Z, m)["energy"]) / (2 * h)
        assert res["forces"][a, k] == pytest.approx(fd, abs=2e-3)
    finer = run(make(method="wb97m-v", basis="def2-svp", nlc_grid_level=2), WATER_Z, WATER_XYZ)["energy"]
    assert abs(finer - res["energy"]) / HARTREE < 1e-6          # the default VV10 grid is converged


def test_wb97m_v_hessian_is_numerical_and_consistent():
    eng = make(method="wb97m-v", basis="def2-svp", density_fit=True)
    assert eng.meta_gga and eng.nonlocal_term
    h2 = [[0, 0, 0], [0, 0, 0.74]]
    hess = run(eng, [1, 1], h2, ["hessian"])["hessian"]
    assert hess.shape == (6, 6) and np.allclose(hess, hess.T, atol=1e-6)
    assert hess[2, 2] > 10.0 and hess[5, 5] > 10.0              # the H-H stretch, a stiff bond (eV/A^2)


def test_other_functionals():
    cam = make(method="other", method_other="  CAM-B3LYP ", basis="def2-svp")        # case and spaces do not matter
    assert cam.xc == "cam-b3lyp" and not cam.meta_gga
    e_cam = run(cam, WATER_Z, WATER_XYZ)["energy"]
    e_b3 = run(make(method="b3lyp", basis="def2-svp"), WATER_Z, WATER_XYZ)["energy"]
    assert abs(e_cam - e_b3) / HARTREE > 1e-3 and -76.5 < e_cam / HARTREE < -76.2     # another functional
    same = make(method="other", method_other="b3lyp", basis="def2-svp")
    assert run(same, WATER_Z, WATER_XYZ)["energy"] == pytest.approx(e_b3, abs=1e-9)
    hf = make(method="other", method_other="hf", basis="sto-3g", density_fit=False)    # hf also works as typed
    assert run(hf, [1, 1], [[0, 0, 0], [0, 0, 0.74]])["energy"] / HARTREE == pytest.approx(-1.116759, abs=2e-6)
    scan = make(method="other", method_other="scan", basis="sto-3g")                  # a meta-GGA: the numerical Hessian
    assert scan.meta_gga
    assert run(scan, [1, 1], [[0, 0, 0], [0, 0, 0.74]], ["hessian"])["hessian"].shape == (6, 6)


def test_an_unknown_functional_is_a_clear_error():
    eng = PySCFEngine()
    with pytest.raises(EngineError) as e:
        eng.load(resolve_options(PySCFEngine, {"method": "other", "method_other": "not-a-functional"}, "load"), Context())
    assert e.value.code == "bad_option" and "not-a-functional" in e.value.message
