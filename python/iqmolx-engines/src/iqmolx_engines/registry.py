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
    from .engines.torchani_engine import ALL_ENGINES as TORCHANI
    from .engines.egret_engine import ALL_ENGINES as EGRET
    classes.extend(TORCHANI)
    classes.extend(EGRET)
    return classes


def _warn(text):
    print("iqmolx-engines: " + text, file=sys.stderr)


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
            _warn("cannot load plugin %s: %s" % (ep.name, exc))
            continue
        if not (isinstance(cls, type) and issubclass(cls, Engine)):
            _warn("plugin %s is not an Engine subclass, ignored" % ep.name)
        elif not cls.name:
            _warn("plugin %s has no engine name, ignored" % ep.name)
        else:
            found.append(cls)
    return found


def engines():
    """Dict name -> Engine subclass.  Built-in names win over plugins (with a warning on stderr)."""
    result = {}
    for cls in _entry_points():
        result[cls.name] = cls
    for cls in _builtin():
        if cls.name in result:
            _warn("the plugin engine '%s' is hidden by the built-in engine of that name" % cls.name)
        result[cls.name] = cls
    return result
