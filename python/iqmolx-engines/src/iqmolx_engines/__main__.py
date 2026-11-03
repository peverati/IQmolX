# Copyright (C) 2026 Roberto Peverati
# Developed with the assistance of Claude (Anthropic).
#
# This file is part of IQMolX, a fork of IQmol (http://iqmol.org). It is free
# software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version. It is
# distributed WITHOUT ANY WARRANTY; see <http://www.gnu.org/licenses/>.

"""Command-line entry point: `python -m iqmolx_engines` or `iqmolx-engines`."""

import argparse
import json
import os
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="iqmolx-engines",
        description="Calculator helper for IQMolX. With no options it serves the JSON-lines "
                    "protocol on stdin/stdout.")
    parser.add_argument("--version", action="store_true", help="print versions and exit")
    parser.add_argument("--describe", action="store_true",
                        help="print the engines and their capabilities as JSON and exit")
    args = parser.parse_args(argv)

    try:
        from . import PROTOCOL, __version__, registry
        from .server import Server
    except Exception as exc:
        print(json.dumps({"event": "fatal", "message": "%s: %s" % (type(exc).__name__, exc)}))
        return 2

    if args.version:
        print("iqmolx-engines %s (protocol %d)" % (__version__, PROTOCOL))
        return 0
    if args.describe:
        print(json.dumps({"engines": [c.describe() for c in registry.engines().values()]},
                         indent=2))
        return 0

    # Only protocol messages may reach stdout.  Engine code that prints goes to stderr, and so does
    # native code (C or Fortran libraries write to file descriptor 1 directly): the protocol keeps
    # a private copy of the real stdout and descriptor 1 is pointed at stderr.
    protocol = os.fdopen(os.dup(1), "w", buffering=1, encoding="utf-8")
    os.dup2(2, 1)
    sys.stdout = sys.stderr
    return Server(sys.stdin, protocol).run()


if __name__ == "__main__":
    sys.exit(main())
