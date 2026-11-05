#!/bin/bash
# Copyright (C) 2026 Roberto Peverati.  Part of IQmol, free software under the GNU General
# Public License, version 3 or (at your option) any later version.  Developed with the
# assistance of Claude (Anthropic).
#
# Packages a built IQMolX for Windows (x64) as a zip file that can be unpacked and run: the
# program, the Qt plugins and every DLL it needs, the data files and the xtb program.
# Run it in an MSYS2 UCRT64 shell after building (see .github/workflows/windows-build.yml).
#
# Usage:  ./package-windows.sh [build-dir] [version]
#   build-dir   default: build
#   version     default: from git
#   XTB_ROOT=<xtb install dir>   where xtb.exe and libxtb come from
#   UV_BIN=<path to uv.exe>      uv, bundled to make the Python environment of the Python
#                                engines (default: uv from the PATH; pip install uv gives one)
#   UV_LICENSE_DIR=<folder>      its license files (default: those of the pip package of uv)
#
# Result:  dist/IQMolX-<version>-windows-x64.zip  with  IQMolX/bin/IQMolX.exe  and  IQMolX/share

set -e

TOPDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
BUILD=${1:-build}
VERSION=${2:-$(git -C "$TOPDIR" describe --tags --match "v[0-9]*" --always 2>/dev/null || echo dev)}
XTB_ROOT=${XTB_ROOT:-$TOPDIR/../xtb-install}
MINGW=${MSYSTEM_PREFIX:-/ucrt64}

EXE="$TOPDIR/$BUILD/bin/IQmol.exe"
[ -f "$EXE" ] || EXE=$(find "$TOPDIR/$BUILD" -name "IQmol.exe" | head -1)
[ -f "$EXE" ] || { echo "IQmol.exe not found in $BUILD"; exit 1; }

NAME="IQMolX-$VERSION-windows-x64"
OUT="$TOPDIR/dist/$NAME/IQMolX"
rm -rf "$TOPDIR/dist/$NAME" "$TOPDIR/dist/$NAME.zip"
mkdir -p "$OUT/bin" "$OUT/share"

cp "$EXE" "$OUT/bin/IQMolX.exe"

# Data files the program looks for next to bin (see Preferences.C)
cp -R "$TOPDIR/share/." "$OUT/share/"
mkdir -p "$OUT/share/openbabel"
cp -R "$TOPDIR/modules/openbabel/data/." "$OUT/share/openbabel/"

# The xtb program and library
if [ -f "$XTB_ROOT/bin/xtb.exe" ]; then
   cp "$XTB_ROOT/bin/xtb.exe" "$OUT/bin/"
   cp "$XTB_ROOT"/bin/*.dll "$OUT/bin/" 2>/dev/null || true
   cp "$XTB_ROOT"/lib/*.dll "$OUT/bin/" 2>/dev/null || true
else
   echo "Warning: no xtb.exe in $XTB_ROOT/bin, it will not be part of the package"
fi

# The Python engines (Calculation > Python Engines): the helper as a wheel, and uv, which makes
# the Python environment for it when the user asks for that.  Weights and PyTorch are never part
# of the package.  Building the wheel needs the internet (setuptools from PyPI).
UV_BIN=${UV_BIN:-$(command -v uv || command -v uv.exe || true)}
if [ -n "$UV_BIN" ] && [ -f "$UV_BIN" ]; then
   mkdir -p "$OUT/python" "$OUT/licenses/uv"
   cp "$UV_BIN" "$OUT/bin/uv.exe"

   # Built from a copy: building writes build/ and *.egg-info next to the sources
   HELPER=$(mktemp -d)/helper-src
   mkdir -p "$HELPER"
   (cd "$TOPDIR/python/iqmolx-engines" && tar --exclude=build --exclude='*.egg-info' --exclude=__pycache__ \
      --exclude=.pytest_cache --exclude=tests -cf - .) | (cd "$HELPER" && tar -xf -)
   "$UV_BIN" build --wheel --out-dir "$OUT/python" "$(cygpath -w "$HELPER" 2>/dev/null || echo "$HELPER")"
   rm -f "$OUT/python/.gitignore"
   ls "$OUT"/python/iqmolx_engines-*.whl > /dev/null

   # uv is MIT or Apache-2.0 licensed: its license texts travel with it
   UV_LICENSES=${UV_LICENSE_DIR:-}
   if [ -z "$UV_LICENSES" ]; then
      for d in "$(dirname "$UV_BIN")"/../Lib/site-packages/uv-*.dist-info/licenses \
               "$(dirname "$UV_BIN")"/../lib/python*/site-packages/uv-*.dist-info/licenses; do
         [ -d "$d" ] && { UV_LICENSES=$d; break; }
      done
   fi
   if [ -z "$UV_LICENSES" ] || [ ! -d "$UV_LICENSES" ]; then
      echo "uv is bundled but its license files were not found; set UV_LICENSE_DIR to a folder with them"
      exit 1
   fi
   cp -R "$UV_LICENSES"/. "$OUT/licenses/uv/"
   echo "uv $("$UV_BIN" --version | awk '{print $2}') from https://github.com/astral-sh/uv" > "$OUT/licenses/uv/SOURCE.txt"
   cp "$TOPDIR/LICENSE" "$OUT/licenses/iqmolx-engines-LICENSE"
else
   echo "Warning: no uv found (set UV_BIN), the Python engines cannot set up their environment"
fi

# Qt plugins that a Widgets program needs, and a qt.conf that points at them
PLUGINS="$MINGW/share/qt5/plugins"
[ -d "$PLUGINS" ] || PLUGINS="$MINGW/lib/qt5/plugins"
mkdir -p "$OUT/bin/plugins"
for d in platforms imageformats sqldrivers printsupport styles iconengines; do
   if [ -d "$PLUGINS/$d" ]; then
      mkdir -p "$OUT/bin/plugins/$d"
      cp "$PLUGINS/$d/"*.dll "$OUT/bin/plugins/$d/" 2>/dev/null || true
   fi
done
rm -f "$OUT"/bin/plugins/*/*d.dll 2>/dev/null || true      # debug builds of the plugins
printf '[Paths]\nPlugins = plugins\n' > "$OUT/bin/qt.conf"

# Every DLL from the MinGW tree that the program, xtb and the plugins need (repeat until nothing
# new is found, since the DLLs have DLLs of their own)
for round in 1 2 3 4 5 6; do
   NEW=0
   for f in "$OUT"/bin/*.exe "$OUT"/bin/*.dll "$OUT"/bin/plugins/*/*.dll; do
      [ -f "$f" ] || continue
      for dll in $(ldd "$f" 2>/dev/null | awk '/=>/ {print $3}' | grep -i "^$MINGW" || true); do
         if [ ! -f "$OUT/bin/$(basename "$dll")" ]; then cp "$dll" "$OUT/bin/"; NEW=1; fi
      done
   done
   [ $NEW -eq 0 ] && break
done

echo "Files in the package: $(find "$OUT" -type f | wc -l), $(du -sh "$OUT" | cut -f1)"

# What would still be loaded from outside the folder?
OUTSIDE=$(for f in "$OUT"/bin/*.exe "$OUT"/bin/*.dll; do ldd "$f" 2>/dev/null | awk '/=>/ {print $3}'; done | grep -i "^$MINGW" | sort -u || true)
if [ -n "$OUTSIDE" ]; then echo "Still linked to libraries outside the package:"; echo "$OUTSIDE"; exit 1; fi

cd "$TOPDIR/dist/$NAME"
zip -qr "../$NAME.zip" IQMolX
echo "Built: $TOPDIR/dist/$NAME.zip ($(du -h "$TOPDIR/dist/$NAME.zip" | cut -f1))"
