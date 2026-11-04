# Copyright (C) 2026 Roberto Peverati
# Developed with the assistance of Claude (Anthropic).
#
# This file is part of IQMolX, a fork of IQmol (http://iqmol.org). It is free
# software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version. It is
# distributed WITHOUT ANY WARRANTY; see <http://www.gnu.org/licenses/>.

"""The EMT plugin of the tutorial (python/examples/iqmolx-engine-emt): the engine, and a real
install that the helper finds and runs.  Needs ASE; skipped without it."""

import json
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest

pytest.importorskip("ase")

HERE = os.path.dirname(__file__)
SRC = os.path.abspath(os.path.join(HERE, "..", "src"))
EXAMPLE = os.path.abspath(os.path.join(HERE, "..", "..", "examples", "iqmolx-engine-emt"))
sys.path.insert(0, os.path.join(EXAMPLE, "src"))

from iqmolx_engine_emt import EmtEngine  # noqa: E402
from iqmolx_engines.base import Context, EngineError, resolve_options  # noqa: E402

CU4 = ([29] * 4, [[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [1.25, 2.165, 0.0], [1.25, 0.722, 2.041]])


def make(**opts):
    eng = EmtEngine()
    eng.load(resolve_options(EmtEngine, opts, "load"), Context())
    return eng


def run(eng, z, xyz, want=("energy", "forces"), **opts):
    return eng.compute(z, np.asarray(xyz, float), 0, 1, list(want),
                       resolve_options(EmtEngine, opts, "compute"), Context())


def test_the_engine_describes_itself():
    d = EmtEngine.describe()
    assert d["name"] == "emt" and d["available"] is True
    assert d["elements"] == [1, 6, 7, 8, 13, 28, 29, 46, 47, 78, 79]
    assert d["charge"] is False and d["multiplicity"] is False
    assert [(o["name"], o["scope"]) for o in d["options"]] == [("asap_cutoff", "load"), ("min_distance", "compute")]
    assert d["requirements"] == ["ase"]


def test_energy_and_forces_of_a_copper_dimer():
    res = run(make(), [29, 29], [[0, 0, 0], [0, 0, 2.4]])
    assert res["energy"] == pytest.approx(3.50440546, abs=1e-6)         # eV, from ASE itself
    assert res["forces"][0] == pytest.approx([0.0, 0.0, 2.35718828], abs=1e-6)
    assert res["forces"][0] == pytest.approx(-res["forces"][1])


def test_forces_are_minus_the_gradient_of_the_energy():
    eng, (z, xyz) = make(), CU4
    f = run(eng, z, xyz)["forces"]
    xyz = np.array(xyz)
    h = 1e-5
    for a in range(4):
        for k in range(3):
            p, m = xyz.copy(), xyz.copy()
            p[a, k] += h
            m[a, k] -= h
            fd = -(run(eng, z, p, ["energy"])["energy"] - run(eng, z, m, ["energy"])["energy"]) / (2 * h)
            assert f[a, k] == pytest.approx(fd, abs=1e-6)


def test_atoms_that_almost_touch_are_refused_and_the_check_can_be_changed():
    eng = make()
    close = [[0, 0, 0], [0, 0, 0.05]]
    with pytest.raises(EngineError) as e:
        run(eng, [29, 29], close)
    assert e.value.code == "engine_error" and "0.05" in e.value.message
    assert run(eng, [29, 29], close, min_distance=0.0)["energy"] == pytest.approx(1263.9, abs=0.1)


def test_the_load_option_reaches_the_calculator():
    assert make(asap_cutoff=True).calc.parameters["asap_cutoff"] is True
    assert make().calc.parameters["asap_cutoff"] is False


def test_what_the_plugin_declares_is_enforced_by_the_helper():
    with pytest.raises(EngineError) as e:
        EmtEngine.validate_system([26, 26], 0, 1, ["energy"])           # iron
    assert e.value.code == "unsupported_element"
    with pytest.raises(EngineError) as e:
        EmtEngine.validate_system([29, 29], -1, 1, ["energy"])
    assert e.value.code == "unsupported_charge"


# --- a real install: the helper finds the plugin in its environment ----------------------------

@pytest.fixture(scope="module")
def installed(tmp_path_factory):
    work = tmp_path_factory.mktemp("emt")
    source = work / "source"
    shutil.copytree(EXAMPLE, source)             # building writes build/ next to the sources
    target = work / "site"
    done = subprocess.run([sys.executable, "-m", "pip", "install", "--quiet", "--no-deps",
                           "--no-build-isolation", "--target", str(target), str(source)],
                          capture_output=True, text=True)
    if done.returncode != 0:
        pytest.skip("cannot pip install the tutorial plugin here: " + done.stderr.strip()[-300:])
    return str(target)


def test_the_helper_finds_the_installed_plugin_and_uses_it(installed):
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([SRC, installed]))
    proc = subprocess.Popen([sys.executable, "-m", "iqmolx_engines"], stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    json.loads(proc.stdout.readline())                       # ready

    def call(i, t, **k):
        proc.stdin.write(json.dumps(dict(id=i, type=t, **k)) + "\n")
        proc.stdin.flush()
        return json.loads(proc.stdout.readline())

    engines = {e["name"]: e for e in call(1, "list_engines")["result"]["engines"]}
    assert engines["emt"]["title"] == "EMT metals (ASE, example plugin)" and engines["emt"]["available"]
    assert call(2, "load", engine="emt", options={"asap_cutoff": False})["ok"]
    res = call(3, "compute", numbers=[29, 29], coordinates=[[0, 0, 0], [0, 0, 2.4]])["result"]
    assert res["energy"] == pytest.approx(3.50440546, abs=1e-6)
    too_close = call(4, "compute", numbers=[29, 29], coordinates=[[0, 0, 0], [0, 0, 0.05]])
    assert too_close["error"]["code"] == "engine_error"
    near_ok = call(5, "compute", numbers=[29, 29], coordinates=[[0, 0, 0], [0, 0, 0.05]],
                   options={"min_distance": 0.0})
    assert near_ok["ok"]
    iron = call(6, "compute", numbers=[26, 26], coordinates=[[0, 0, 0], [0, 0, 2.0]])
    assert iron["error"]["code"] == "unsupported_element"
    proc.stdin.close()
    assert proc.wait(timeout=10) == 0
