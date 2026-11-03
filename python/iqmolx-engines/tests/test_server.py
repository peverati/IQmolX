import json
import os
import subprocess
import sys
import time

import pytest

SRC = os.path.join(os.path.dirname(__file__), "..", "src")
WATER = {"numbers": [8, 1, 1],
         "coordinates": [[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]]}


class Helper:
    """The helper as IQMolX will run it: a child process speaking JSON lines."""

    def __init__(self):
        env = dict(os.environ, PYTHONPATH=SRC)
        self.proc = subprocess.Popen([sys.executable, "-m", "iqmolx_engines"],
                                     stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE, text=True, env=env)
        self.ready = self.read()

    def send(self, **msg):
        self.proc.stdin.write(json.dumps(msg) + "\n")
        self.proc.stdin.flush()

    def send_raw(self, text):
        self.proc.stdin.write(text + "\n")
        self.proc.stdin.flush()

    def read(self):
        line = self.proc.stdout.readline()
        assert line, "helper closed stdout: " + self.proc.stderr.read()
        return json.loads(line)

    def final(self, rid):
        """Read until the final message of request `rid`; return (events, message)."""
        events = []
        while True:
            msg = self.read()
            if msg.get("id") == rid and "event" not in msg:
                return events, msg
            events.append(msg)

    def call(self, rid, type, **fields):
        self.send(id=rid, type=type, **fields)
        return self.final(rid)[1]

    def close(self):
        self.proc.stdin.close()
        assert self.proc.wait(timeout=10) == 0


@pytest.fixture
def helper():
    h = Helper()
    yield h
    if h.proc.poll() is None:
        h.close()


def test_ready_and_hello(helper):
    assert helper.ready["event"] == "ready" and helper.ready["protocol"] == 1
    r = helper.call(1, "hello", protocol=1, client="test")
    assert r["ok"] and r["result"]["protocol"] == 1
    r = helper.call(2, "hello", protocol=7)
    assert r["result"]["protocol"] == 1
    r = helper.call(3, "hello", protocol=0)
    assert r["error"]["code"] == "unsupported_protocol"


def test_list_engines_has_capabilities(helper):
    r = helper.call(1, "list_engines")
    demo = next(e for e in r["result"]["engines"] if e["name"] == "demo")
    assert demo["available"] and "hessian" in demo["properties"]
    assert {o["name"] for o in demo["options"]} == {"well_depth", "alpha", "delay_ms"}
    assert demo["requirements"] == []
    by_name = {e["name"]: e for e in r["result"]["engines"]}
    assert by_name["torchani"]["requirements"] == ["torch", "torchani"]
    assert by_name["ase"]["requirements"] == ["ase"]


def test_compute_flow(helper):
    assert helper.call(1, "compute", **WATER)["error"]["code"] == "not_loaded"
    r = helper.call(2, "load", engine="demo", options={"well_depth": 2.0})
    assert r["ok"] and r["result"]["engine"] == "demo"
    r = helper.call(3, "compute", want=["energy", "forces", "hessian"], **WATER)["result"]
    assert isinstance(r["energy"], float)
    assert len(r["forces"]) == 3 and len(r["hessian"]) == 9
    assert set(r) == {"energy", "forces", "hessian"}
    r = helper.call(4, "shutdown")
    assert r["ok"]
    assert helper.proc.wait(timeout=10) == 0


def test_errors_keep_helper_alive(helper):
    helper.send_raw("this is not json")
    assert helper.read()["error"]["code"] == "bad_request"
    assert helper.call(1, "nonsense")["error"]["code"] == "bad_request"
    assert helper.call(2, "load", engine="nope")["error"]["code"] == "unknown_engine"
    assert helper.call(3, "load", engine="demo", options={"alpha": "x"}
                       )["error"]["code"] == "bad_option"
    helper.call(4, "load", engine="demo")
    bad = dict(WATER, numbers=[92, 1, 1])
    assert helper.call(5, "compute", **bad)["error"]["code"] == "unsupported_element"
    assert helper.call(6, "compute", want=["spin"], **WATER)["error"]["code"] \
        == "unsupported_property"
    assert helper.call(7, "compute", numbers=[1], coordinates=[[0, 0]]
                       )["error"]["code"] == "bad_request"
    assert helper.call(8, "compute", **WATER)["ok"]


def test_cancel_running_computation(helper):
    helper.call(1, "load", engine="demo")
    t0 = time.monotonic()
    helper.send(id=2, type="compute", options={"delay_ms": 30000}, **WATER)
    time.sleep(0.2)
    r = helper.call(3, "cancel", target=2)
    assert r["result"]["found"] is True
    _, final = helper.final(2)
    assert final["error"]["code"] == "cancelled"
    assert time.monotonic() - t0 < 10
    assert helper.call(4, "cancel", target=2)["result"]["found"] is False
    assert helper.call(5, "compute", **WATER)["ok"]


def test_cancel_queued_request_never_runs(helper):
    helper.call(1, "load", engine="demo")
    helper.send(id=2, type="compute", options={"delay_ms": 600}, **WATER)
    helper.send(id=3, type="compute", options={"delay_ms": 30000}, **WATER)
    helper.send(id=4, type="cancel", target=3)
    seen = {}
    while len(seen) < 3:
        msg = helper.read()
        seen[msg["id"]] = msg
    assert seen[4]["result"]["found"] is True
    assert seen[3]["error"]["code"] == "cancelled"
    assert seen[2]["ok"]


def test_stdout_stays_clean_when_engine_prints(helper, tmp_path):
    # An engine that prints must not corrupt the protocol stream.
    plugin = tmp_path / "noisy.py"
    plugin.write_text(
        "from iqmolx_engines.engines.demo import DemoEngine\n"
        "class Noisy(DemoEngine):\n"
        "    name = 'noisy'\n"
        "    def compute(self, *a):\n"
        "        print('hello from the engine')\n"
        "        return super().compute(*a)\n")
    helper.close()
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([SRC, str(tmp_path)]))
    code = ("import sys, io\nsys.argv=['x']\n"
            "import noisy\n"
            "from iqmolx_engines import registry\n"
            "registry._builtin = lambda c=registry._builtin: c() + [noisy.Noisy]\n"
            "from iqmolx_engines.__main__ import main\nsys.exit(main())\n")
    proc = subprocess.Popen([sys.executable, "-c", code], stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                            env=env)
    msgs = [{"id": 1, "type": "load", "engine": "noisy"},
            dict(id=2, type="compute", **WATER)]
    out, err = proc.communicate("".join(json.dumps(m) + "\n" for m in msgs), timeout=20)
    lines = [json.loads(line) for line in out.splitlines()]
    assert lines[0]["event"] == "ready" and lines[2]["ok"]
    assert "hello from the engine" in err


def test_cli_describe_and_version():
    env = dict(os.environ, PYTHONPATH=SRC)
    out = subprocess.run([sys.executable, "-m", "iqmolx_engines", "--describe"],
                         capture_output=True, text=True, env=env, check=True).stdout
    assert any(e["name"] == "demo" for e in json.loads(out)["engines"])
    out = subprocess.run([sys.executable, "-m", "iqmolx_engines", "--version"],
                         capture_output=True, text=True, env=env, check=True).stdout
    assert "protocol 1" in out


def test_required_option_is_listed(helper):
    engines = {e["name"]: e for e in helper.call(1, "list_engines")["result"]["engines"]}
    calc = next(o for o in engines["ase"]["options"] if o["name"] == "calculator")
    assert calc["required"] is True and calc["placeholder"] == "ase.calculators.emt:EMT"


def test_empty_required_option_is_refused_by_the_helper(helper):
    engines = {e["name"]: e for e in helper.call(1, "list_engines")["result"]["engines"]}
    if not engines["ase"]["available"]:
        pytest.skip("needs ASE: the helper answers engine_unavailable first")
    for options in ({}, {"calculator": ""}, {"calculator": "   "}):
        err = helper.call(2, "load", engine="ase", options=options)["error"]
        assert err["code"] == "bad_option" and "Calculator class" in err["message"]
