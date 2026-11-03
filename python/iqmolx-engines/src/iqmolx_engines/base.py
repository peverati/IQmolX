# Copyright (C) 2026 Roberto Peverati
# Developed with the assistance of Claude (Anthropic).
#
# This file is part of IQMolX, a fork of IQmol (http://iqmol.org). It is free
# software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version. It is
# distributed WITHOUT ANY WARRANTY; see <http://www.gnu.org/licenses/>.

"""Engine base class and the errors that map to protocol error codes."""

import numpy as np


class EngineError(Exception):
    """An error with a protocol error code (see PROTOCOL.md)."""

    def __init__(self, code, message, detail=""):
        super().__init__(message)
        self.code = code
        self.message = message
        self.detail = detail


class Cancelled(EngineError):
    def __init__(self, message="Cancelled"):
        super().__init__("cancelled", message)


class Context:
    """Passed to engines so they can report progress and notice cancellation."""

    def __init__(self, progress=None, cancelled=None):
        self._progress = progress
        self._cancelled = cancelled

    def progress(self, fraction=None, message=""):
        if self._progress:
            self._progress(fraction, message)

    def cancelled(self):
        return bool(self._cancelled and self._cancelled())

    def check_cancel(self):
        if self.cancelled():
            raise Cancelled()


def option(name, title, type, default, scope="load", **extra):
    """Build an option description for the capabilities.  `extra` may hold min, max, choices,
    description, `required=True` (the option must not be empty: the window keeps Run off until it
    is filled in) and `placeholder` (a hint shown in an empty text field)."""
    d = {"name": name, "title": title, "type": type, "default": default, "scope": scope}
    d.update(extra)
    return d


class Engine:
    """Base class of all engines. Subclasses set the class attributes and implement load()
    and compute(). Everything in `describe()` is what IQMolX sees."""

    name = ""
    title = ""
    description = ""
    version = "0.1.0"
    license = ""
    citation = ""
    elements = None            # list of atomic numbers, None for all
    charge = True
    multiplicity = True
    properties = ("energy", "forces")
    options = ()               # tuple of option() dicts
    requires_download = False
    requires_login = False
    requirements = ()          # pip requirement strings the engine needs, for the setup step

    @classmethod
    def unavailable_reason(cls):
        """Empty string when the engine can run; otherwise what is missing."""
        return ""

    @classmethod
    def describe(cls):
        reason = cls.unavailable_reason()
        return {
            "name": cls.name, "title": cls.title, "description": cls.description,
            "version": cls.version, "license": cls.license, "citation": cls.citation,
            "available": not reason, "unavailable_reason": reason,
            "elements": None if cls.elements is None else list(cls.elements),
            "charge": cls.charge, "multiplicity": cls.multiplicity,
            "properties": list(cls.properties),
            "options": [dict(o) for o in cls.options],
            "requires_download": cls.requires_download,
            "requires_login": cls.requires_login,
            "requirements": list(cls.requirements),
        }

    def load(self, options, ctx):
        """Load the model. `options` is complete (defaults filled in). Return an info dict."""
        return {}

    def compute(self, numbers, coordinates, charge, multiplicity, want, options, ctx):
        """Return a dict with the wanted properties as numpy arrays or floats."""
        raise NotImplementedError

    # -- validation shared by the server ------------------------------------------------

    @classmethod
    def validate_system(cls, numbers, charge, multiplicity, want):
        if cls.elements is not None:
            bad = sorted({int(z) for z in numbers} - set(cls.elements))
            if bad:
                raise EngineError("unsupported_element",
                                  "Unsupported element(s): %s" % ", ".join(map(str, bad)),
                                  detail="supported: %s" % list(cls.elements))
        if charge != 0 and not cls.charge:
            raise EngineError("unsupported_charge",
                              "%s supports neutral molecules only" % cls.title)
        if multiplicity != 1 and not cls.multiplicity:
            raise EngineError("unsupported_multiplicity",
                              "%s supports closed-shell (multiplicity 1) only" % cls.title)
        missing = [w for w in want if w not in cls.properties]
        if missing:
            raise EngineError("unsupported_property",
                              "Not available from %s: %s" % (cls.title, ", ".join(missing)))


def check_option_value(opt, value):
    """Return the value coerced to the option's type, or raise bad_option."""
    t, name = opt["type"], opt["name"]
    try:
        if t == "bool":
            if not isinstance(value, bool):
                raise ValueError
        elif t == "int":
            if isinstance(value, bool) or not isinstance(value, int):
                raise ValueError
        elif t == "float":
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError
            value = float(value)
        elif t == "string":
            if not isinstance(value, str):
                raise ValueError
        elif t == "choice":
            if value not in opt.get("choices", []):
                raise ValueError
        else:
            raise ValueError
    except ValueError:
        raise EngineError("bad_option", "Option '%s' must be of type %s" % (name, t)) from None
    if t in ("int", "float"):
        if "min" in opt and value < opt["min"] or "max" in opt and value > opt["max"]:
            raise EngineError("bad_option", "Option '%s' is out of range" % name)
    return value


def resolve_options(engine_cls, given, scope):
    """Fill defaults and validate the options of one scope ('load' or 'compute')."""
    known = {o["name"]: o for o in engine_cls.options if o.get("scope", "load") == scope}
    given = given or {}
    if not isinstance(given, dict):
        raise EngineError("bad_request", "'options' must be an object")
    for key in given:
        if key not in known:
            raise EngineError("bad_option", "Unknown %s option '%s'" % (scope, key))
    values = {n: check_option_value(o, given[n]) if n in given else o["default"]
              for n, o in known.items()}
    for n, o in known.items():
        value = values[n]
        if o.get("required") and (value is None or (isinstance(value, str) and not value.strip())):
            raise EngineError("bad_option", "'%s' is required" % o.get("title", n))
    return values


def as_system(numbers, coordinates):
    """Validate and convert the molecule of a compute request to numpy arrays."""
    try:
        z = np.asarray(numbers, dtype=int)
        xyz = np.asarray(coordinates, dtype=float)
    except (TypeError, ValueError):
        raise EngineError("bad_request", "'numbers' and 'coordinates' must be numeric") from None
    if z.ndim != 1 or len(z) == 0 or xyz.shape != (len(z), 3):
        raise EngineError("bad_request", "Expected N atomic numbers and N x 3 coordinates")
    if not np.all(np.isfinite(xyz)):
        raise EngineError("bad_request", "Coordinates must be finite")
    if np.any(z < 1):
        raise EngineError("bad_request", "Atomic numbers must be positive")
    return z, xyz
