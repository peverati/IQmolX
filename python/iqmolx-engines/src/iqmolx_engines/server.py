# Copyright (C) 2026 Roberto Peverati
# Developed with the assistance of Claude (Anthropic).
#
# This file is part of IQMolX, a fork of IQmol (http://iqmol.org). It is free
# software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version. It is
# distributed WITHOUT ANY WARRANTY; see <http://www.gnu.org/licenses/>.

"""The protocol loop: reads JSON lines, runs requests one at a time, answers each once.

A reader thread parses requests and handles `cancel` immediately, so a running computation
can be flagged while the worker (the calling thread) is busy. See PROTOCOL.md."""

import json
import platform
import queue
import sys
import threading
import traceback

import numpy as np

from . import PROTOCOL, __version__
from .base import (Cancelled, Context, EngineError, as_system, resolve_options)
from . import registry


def _jsonable(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    raise TypeError("not JSON serializable: %r" % type(obj))


class Server:
    def __init__(self, instream, outstream):
        self.inp = instream
        self.out = outstream
        self.write_lock = threading.Lock()
        self.queue = queue.Queue()
        self.lock = threading.Lock()
        self.pending = set()       # ids queued or running
        self.cancelled = set()     # ids flagged by cancel
        self.engine = None
        self.engine_cls = None
        self.running = True

    # -- output -------------------------------------------------------------------------

    def send(self, msg):
        line = json.dumps(msg, default=_jsonable, allow_nan=False, separators=(",", ":"))
        with self.write_lock:
            self.out.write(line + "\n")
            self.out.flush()

    def reply(self, rid, result):
        self.send({"id": rid, "ok": True, "result": result})

    def fail(self, rid, code, message, detail=""):
        err = {"code": code, "message": message}
        if detail:
            err["detail"] = detail
        self.send({"id": rid, "ok": False, "error": err})

    # -- reader thread ------------------------------------------------------------------

    def read_loop(self):
        try:
            for line in self.inp:
                line = line.strip()
                if line:
                    self.accept(line)
        finally:
            self.queue.put(None)

    def accept(self, line):
        try:
            req = json.loads(line)
        except ValueError:
            return self.fail(None, "bad_request", "Invalid JSON")
        if not isinstance(req, dict):
            return self.fail(None, "bad_request", "A request must be a JSON object")
        rid = req.get("id")
        if isinstance(rid, bool) or not isinstance(rid, (int, str)):
            return self.fail(None, "bad_request", "'id' must be an integer or a string")
        if not isinstance(req.get("type"), str):
            return self.fail(rid, "bad_request", "'type' must be a string")
        with self.lock:
            if rid in self.pending:
                return self.fail(rid, "bad_request", "Duplicate id %r" % (rid,))
            if req["type"] == "cancel":
                target = req.get("target")
                found = isinstance(target, (int, str)) and target in self.pending
                if found:
                    self.cancelled.add(target)
                return self.reply(rid, {"found": found})
            self.pending.add(rid)
        self.queue.put(req)

    # -- worker -------------------------------------------------------------------------

    def run(self):
        self.send({"event": "ready", "protocol": PROTOCOL, "helper": __version__,
                   "python": platform.python_version()})
        threading.Thread(target=self.read_loop, daemon=True).start()
        while self.running:
            req = self.queue.get()
            if req is None:
                break
            self.process(req)
        return 0

    def is_cancelled(self, rid):
        with self.lock:
            return rid in self.cancelled

    def process(self, req):
        rid = req["id"]
        try:
            if self.is_cancelled(rid):
                raise Cancelled()
            ctx = Context(progress=lambda f, m: self.progress(rid, f, m),
                          cancelled=lambda: self.is_cancelled(rid))
            handler = getattr(self, "do_" + req["type"], None)
            if handler is None:
                raise EngineError("bad_request", "Unknown request type '%s'" % req["type"])
            result = handler(req, ctx)
            if self.is_cancelled(rid):
                raise Cancelled()
            self.reply(rid, result)
        except EngineError as exc:
            self.fail(rid, exc.code, exc.message, exc.detail)
        except Exception as exc:
            self.fail(rid, "internal", "%s: %s" % (type(exc).__name__, exc),
                      traceback.format_exc())
        finally:
            with self.lock:
                self.pending.discard(rid)
                self.cancelled.discard(rid)

    def progress(self, rid, fraction, message):
        msg = {"id": rid, "event": "progress", "message": message}
        if fraction is not None:
            msg["fraction"] = float(fraction)
        self.send(msg)

    # -- handlers -----------------------------------------------------------------------

    def do_hello(self, req, ctx):
        wanted = req.get("protocol")
        if isinstance(wanted, bool) or not isinstance(wanted, int) or wanted < 1:
            raise EngineError("unsupported_protocol", "No common protocol version",
                              "helper speaks up to %d" % PROTOCOL)
        return {"protocol": min(wanted, PROTOCOL), "helper": __version__,
                "python": platform.python_version(),
                "platform": "%s-%s" % (sys.platform, platform.machine())}

    def do_list_engines(self, req, ctx):
        return {"engines": [cls.describe() for cls in registry.engines().values()]}

    def do_load(self, req, ctx):
        name = req.get("engine")
        cls = registry.engines().get(name)
        if cls is None:
            raise EngineError("unknown_engine", "Unknown engine '%s'" % (name,))
        reason = cls.unavailable_reason()
        if reason:
            raise EngineError("engine_unavailable", reason)
        options = resolve_options(cls, req.get("options"), "load")
        engine = cls()
        info = self.call_engine(engine.load, options, ctx)
        self.engine, self.engine_cls = engine, cls
        return {"engine": name, "info": info or {}}

    def do_compute(self, req, ctx):
        if self.engine is None:
            raise EngineError("not_loaded", "No engine is loaded")
        cls = self.engine_cls
        numbers, xyz = as_system(req.get("numbers"), req.get("coordinates"))
        charge, mult = req.get("charge", 0), req.get("multiplicity", 1)
        for value, what in ((charge, "charge"), (mult, "multiplicity")):
            if isinstance(value, bool) or not isinstance(value, int):
                raise EngineError("bad_request", "'%s' must be an integer" % what)
        if mult < 1:
            raise EngineError("bad_request", "'multiplicity' must be at least 1")
        want = req.get("want", ["energy", "forces"])
        if not isinstance(want, list) or not all(isinstance(w, str) for w in want):
            raise EngineError("bad_request", "'want' must be a list of strings")
        cls.validate_system(numbers, charge, mult, want)
        options = resolve_options(cls, req.get("options"), "compute")
        raw = self.call_engine(self.engine.compute, numbers, xyz, charge, mult, want,
                               options, ctx)
        return self.check_result(raw, want, len(numbers))

    def do_shutdown(self, req, ctx):
        self.running = False
        return {}

    # -- helpers ------------------------------------------------------------------------

    @staticmethod
    def call_engine(func, *args):
        """Run engine code; its exceptions become engine_error, the helper stays alive."""
        try:
            return func(*args)
        except EngineError:
            raise
        except Exception as exc:
            raise EngineError("engine_error", "%s: %s" % (type(exc).__name__, exc),
                              traceback.format_exc()) from None

    @staticmethod
    def check_result(raw, want, n):
        shapes = {"energy": (), "forces": (n, 3), "hessian": (3 * n, 3 * n),
                  "dipole": (3,), "charges": (n,)}
        out = {}
        for key in want:
            if key not in raw:
                raise EngineError("engine_error", "Engine did not return '%s'" % key)
            value = np.asarray(raw[key], dtype=float)
            if value.shape != shapes.get(key, value.shape):
                raise EngineError("engine_error", "'%s' has shape %s, expected %s"
                                  % (key, value.shape, shapes[key]))
            if not np.all(np.isfinite(value)):
                raise EngineError("engine_error", "'%s' contains non-finite numbers" % key)
            out[key] = float(value) if key == "energy" else value
        return out
