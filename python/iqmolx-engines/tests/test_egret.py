import os
import re

import numpy as np
import pytest

from iqmolx_engines.base import Context, EngineError, resolve_options
from iqmolx_engines.engines import egret_engine as egret
from iqmolx_engines.engines.egret_engine import ALL_ENGINES

WATER_Z = [8, 1, 1]
WATER_XYZ = np.array([[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]])


def test_every_egret_model_describes_itself_properly():
    names = [c.name for c in ALL_ENGINES]
    assert names == ["egret1", "egret1e", "egret1t", "egret1m", "egret1s"]
    for cls in ALL_ENGINES:
        d = cls.describe()
        assert d["title"] and d["description"] and "MIT" in d["license"] and "Rowan" in d["citation"]
        assert d["charge"] is False and d["multiplicity"] is False
        assert d["requires_download"] is True and d["requirements"] == ["mace-torch>=0.3.14"]
        assert {"energy", "forces", "hessian"} <= set(d["properties"])
        assert re.fullmatch(r"[0-9a-f]{40}", cls.blob_sha1) and cls.model_size > 1_000_000
        assert cls.model_file.endswith(".model")


def test_elements_follow_the_preprint():
    assert egret.Egret1Engine.elements == (1, 6, 7, 8, 9, 15, 16, 17, 35, 53)   # H C N O F P S Cl Br I
    assert egret.Egret1tEngine.elements == egret.Egret1Engine.elements
    assert set(egret.Egret1eEngine.elements) == set(egret.Egret1Engine.elements) | {14}   # adds silicon
    assert 14 in egret.Egret1eEngine.elements and 14 not in egret.Egret1Engine.elements
    for cls, z in ((egret.Egret1Engine, [14, 1]), (egret.Egret1Engine, [11, 17])):    # Si, Na
        with pytest.raises(EngineError) as e:
            cls.validate_system(z, 0, 1, ["energy"])
        assert e.value.code == "unsupported_element"
    egret.Egret1eEngine.validate_system([14, 1, 1, 1, 1], 0, 1, ["energy"])          # SiH4 is fine for 1e


def test_charge_and_spin_are_refused():
    for charge, mult, code in ((1, 1, "unsupported_charge"), (0, 3, "unsupported_multiplicity")):
        with pytest.raises(EngineError) as e:
            egret.Egret1Engine.validate_system(WATER_Z, charge, mult, ["energy"])
        assert e.value.code == code


def test_weights_come_from_a_pinned_commit_and_the_user_cache(monkeypatch, tmp_path):
    assert re.fullmatch(r"[0-9a-f]{40}", egret.COMMIT)
    url = egret.URL % (egret.COMMIT, egret.Egret1Engine.model_file)
    assert url == ("https://raw.githubusercontent.com/rowansci/egret-public/%s/compiled_models/EGRET_1.model"
                   % egret.COMMIT)
    monkeypatch.setenv("IQMOLX_ENGINES_CACHE", str(tmp_path))
    assert egret.Egret1Engine.weights_path() == str(tmp_path / "egret" / egret.COMMIT[:10] / "EGRET_1.model")


def test_unavailable_without_mace():
    reason = egret.Egret1Engine.unavailable_reason()
    try:
        import mace  # noqa: F401
        assert reason == ""
    except ImportError:
        assert "mace-torch" in reason


def _models_asked_for():
    value = os.environ.get("IQMOLX_TEST_MODELS", "")
    if os.environ.get("IQMOLX_TEST_DOWNLOADS") != "1" or egret.Egret1Engine.unavailable_reason():
        return []
    names = [c.name for c in ALL_ENGINES] if value == "all" else [v for v in value.split(",") if v]
    return [c for c in ALL_ENGINES if c.name in names]


@pytest.mark.parametrize("cls", _models_asked_for(), ids=lambda c: c.name)
def test_each_model_gives_sensible_water(cls):
    eng = cls()
    eng.load(resolve_options(cls, {}, "load"), Context())
    res = eng.compute(WATER_Z, WATER_XYZ, 0, 1, ["energy", "forces", "hessian"], {}, Context())
    h = 1e-4
    for a, k in ((0, 2), (1, 1), (2, 1)):
        p, m = WATER_XYZ.copy(), WATER_XYZ.copy()
        p[a, k] += h
        m[a, k] -= h
        ep = eng.compute(WATER_Z, p, 0, 1, ["energy"], {}, Context())["energy"]
        em = eng.compute(WATER_Z, m, 0, 1, ["energy"], {}, Context())["energy"]
        assert res["forces"][a, k] == pytest.approx(-(ep - em) / (2 * h), abs=2e-4)
    assert res["hessian"].shape == (9, 9) and np.allclose(res["hessian"], res["hessian"].T, atol=1e-5)
    assert np.allclose(res["forces"].sum(axis=0), 0.0, atol=2e-3)
    moved = eng.compute(WATER_Z, WATER_XYZ + [2.0, -1.0, 3.0], 0, 1, ["energy"], {}, Context())["energy"]
    assert moved == pytest.approx(res["energy"], abs=1e-6)
