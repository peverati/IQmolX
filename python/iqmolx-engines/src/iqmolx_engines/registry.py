# Copyright (C) 2026 Roberto Peverati
# Developed with the assistance of Claude (Anthropic).
#
# This file is part of IQMolX, a fork of IQmol (http://iqmol.org). It is free
# software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version. It is
# distributed WITHOUT ANY WARRANTY; see <http://www.gnu.org/licenses/>.

"""Engine registry: built-in engines plus those found through the `iqmolx_engines`
entry-point group."""

import sys

from .base import Engine


def _builtin():
    from .engines.demo import DemoEngine
    classes = [DemoEngine]
    try:
        from .engines.ase_adapter import AseEngine
        classes.append(AseEngine)
    except ImportError:
        pass
    return classes


def _entry_points():
    try:
        from importlib.metadata import entry_points
        eps = entry_points()
        group = eps.select(group="iqmolx_engines") if hasattr(eps, "select") \
            else eps.get("iqmolx_engines", [])
    except Exception:
        return []
    found = []
    for ep in group:
        try:
            cls = ep.load()
        except Exception as exc:
            print("iqmolx-engines: cannot load plugin %s: %s" % (ep.name, exc),
                  file=sys.stderr)
            continue
        if isinstance(cls, type) and issubclass(cls, Engine):
            found.append(cls)
    return found


def engines():
    """Dict name -> Engine subclass. Built-in names win over plugins."""
    result = {}
    for cls in _entry_points() + _builtin():
        result[cls.name] = cls
    return result
