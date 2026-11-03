"""The sample plugin of the developer guide (python/examples/iqmolx-engine-lennard-jones): the
engine itself, the plugin discovery of the registry, and a real install that the helper finds."""

import json
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest

HERE = os.path.dirname(__file__)
SRC = os.path.abspath(os.path.join(HERE, "..", "src"))
EXAMPLE = os.path.abspath(os.path.join(HERE, "..", "..", "examples", "iqmolx-engine-lennard-jones"))
sys.path.insert(0, os.path.join(EXAMPLE, "src"))

from iqmolx_engine_lj import LennardJonesEngine, PARAMETERS  # noqa: E402
from iqmolx_engines import registry  # noqa: E402
from iqmolx_engines.base import Context, EngineError, resolve_options  # noqa: E402


def make(**opts):
    eng = LennardJonesEngine()
    eng.load(resolve_options(LennardJonesEngine, opts, "load"), Context())
    return eng


def run(eng, z, xyz, want=("energy", "forces")):
    return eng.compute(z, np.asarray(xyz, float), 0, 1, list(want), {}, Context())


ARGON3 = ([18, 18, 18], [[0.0, 0.0, 0.0], [3.9, 0.0, 0.0], [1.7, 3.3, 0.4]])


def test_argon_dimer_minimum_is_minus_epsilon():
    r0 = 2 ** (1 / 6) * PARAMETERS[18][0]
    res = run(make(), [18, 18], [[0, 0, 0], [0, 0, r0]])
    assert res["energy"] == pytest.approx(-PARAMETERS[18][1])
    assert np.allclose(res["forces"], 0.0, atol=1e-12)


def test_forces_are_minus_the_gradient_of_the_energy():
    eng, (z, xyz) = make(), ARGON3
    f = run(eng, z, xyz)["forces"]
    xyz = np.array(xyz)
    h = 1e-6
    for a in range(3):
        for k in range(3):
            p, m = xyz.copy(), xyz.copy()
            p[a, k] += h
            m[a, k] -= h
            fd = -(run(eng, z, p, ["energy"])["energy"] - run(eng, z, m, ["energy"])["energy"]) / (2 * h)
            assert f[a, k] == pytest.approx(fd, abs=1e-7)
    assert np.allclose(f.sum(axis=0), 0.0, atol=1e-12)


def test_mixed_pair_and_option():
    z, xyz = [2, 54], [[0, 0, 0], [0, 0, 4.0]]
    e1 = run(make(), z, xyz)["energy"]
    e2 = run(make(epsilon_scale=2.0), z, xyz)["energy"]
    assert e2 == pytest.approx(2.0 * e1)


def test_what_the_plugin_declares_is_enforced_by_the_helper():
    with pytest.raises(EngineError) as e:
        LennardJonesEngine.validate_system([8, 1, 1], 0, 1, ["energy"])
    assert e.value.code == "unsupported_element"
    with pytest.raises(EngineError) as e:
        LennardJonesEngine.validate_system([18, 18], 1, 1, ["energy"])
    assert e.value.code == "unsupported_charge"
    with pytest.raises(EngineError) as e:
        LennardJonesEngine.validate_system([18, 18], 0, 1, ["hessian"])
    assert e.value.code == "unsupported_property"


def test_cancel_stops_a_long_computation():
    from iqmolx_engines.base import Cancelled
    eng = make()
    xyz = np.random.default_rng(1).uniform(0, 30, size=(40, 3))
    with pytest.raises(Cancelled):
        eng.compute([18] * 40, xyz, 0, 1, ["energy"], {}, Context(cancelled=lambda: True))


# --- the registry: what it does with plugins, with fake entry points --------------------------

class FakeEntryPoint:
    def __init__(self, name, target):
        self.name, self._target = name, target

    def load(self):
        if isinstance(self._target, Exception):
            raise self._target
        return self._target


class FakeEntryPoints:
    def __init__(self, items):
        self.items = items

    def select(self, group):
        return self.items if group == "iqmolx_engines" else []


def with_plugins(monkeypatch, items):
    import importlib.metadata as metadata
    monkeypatch.setattr(metadata, "entry_points", lambda: FakeEntryPoints(items))


def test_registry_adds_a_plugin(monkeypatch):
    with_plugins(monkeypatch, [FakeEntryPoint("lj", LennardJonesEngine)])
    assert registry.engines()["lj"] is LennardJonesEngine


def test_registry_ignores_broken_plugins_and_says_why(monkeypatch, capsys):
    class NoName(LennardJonesEngine):
        name = ""

    class Impostor(LennardJonesEngine):
        name = "demo"

    with_plugins(monkeypatch, [
        FakeEntryPoint("broken", ImportError("no module named foo")),
        FakeEntryPoint("not_an_engine", object),
        FakeEntryPoint("nameless", NoName),
        FakeEntryPoint("impostor", Impostor),
        FakeEntryPoint("good", LennardJonesEngine),
    ])
    found = registry.engines()
    err = capsys.readouterr().err
    assert "cannot load plugin broken: no module named foo" in err
    assert "plugin not_an_engine is not an Engine subclass" in err
    assert "plugin nameless has no engine name" in err
    assert "plugin engine 'demo' is hidden by the built-in engine" in err
    assert found["lj"] is LennardJonesEngine
    assert found["demo"] is not Impostor          # built-in names win


# --- a real install: the helper finds the plugin in its environment ----------------------------

@pytest.fixture(scope="module")
def installed(tmp_path_factory):
    work = tmp_path_factory.mktemp("plugin")
    source = work / "source"
    shutil.copytree(EXAMPLE, source)             # building writes build/ next to the sources
    target = work / "site"
    done = subprocess.run([sys.executable, "-m", "pip", "install", "--quiet", "--no-deps",
                           "--no-build-isolation", "--target", str(target), str(source)],
                          capture_output=True, text=True)
    if done.returncode != 0:
        pytest.skip("cannot pip install the sample plugin here: " + done.stderr.strip()[-300:])
    return str(target)


def helper(installed):
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([SRC, installed]))
    return subprocess.Popen([sys.executable, "-m", "iqmolx_engines"], stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)


def test_the_helper_finds_the_installed_plugin_and_uses_it(installed):
    proc = helper(installed)
    json.loads(proc.stdout.readline())                       # ready

    def call(i, t, **k):
        proc.stdin.write(json.dumps(dict(id=i, type=t, **k)) + "\n")
        proc.stdin.flush()
        return json.loads(proc.stdout.readline())

    engines = {e["name"]: e for e in call(1, "list_engines")["result"]["engines"]}
    lj = engines["lj"]
    assert lj["title"] == "Lennard-Jones (example plugin)" and lj["available"] is True
    assert lj["elements"] == [2, 10, 18, 36, 54] and lj["charge"] is False
    assert [o["name"] for o in lj["options"]] == ["epsilon_scale"]

    assert call(2, "load", engine="lj")["ok"]
    r0 = 2 ** (1 / 6) * 3.40
    res = call(3, "compute", numbers=[18, 18], coordinates=[[0, 0, 0], [0, 0, r0]])["result"]
    assert res["energy"] == pytest.approx(-0.01032)
    bad = call(4, "compute", numbers=[8, 1, 1], coordinates=[[0, 0, 0], [0, 0, 1], [0, 1, 0]])
    assert bad["error"]["code"] == "unsupported_element"
    proc.stdin.close()
    assert proc.wait(timeout=10) == 0
