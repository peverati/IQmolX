#!/bin/bash
# Copyright (C) 2026 Roberto Peverati.  Part of IQmol, free software under the GNU General
# Public License, version 3 or (at your option) any later version.  Developed with the
# assistance of Claude (Anthropic).
#
# Packages a built IQMolX.app (Apple Silicon) so that it runs on another Mac: the Qt
# frameworks and every other library it needs are copied into the app, the xtb program is
# added if it can be found, the app is signed ad hoc and put into a disk image.
#
# Usage:  ./package-macos.sh [build-dir] [version]
#   build-dir   default: build-clean (the folder build-macos.sh made)
#   version     default: from git
#   XTB_ROOT=<xtb install dir>   where xtb and libxtb come from (the one used for the build)
#
# The result is   dist/IQMolX-<version>-macOS-arm64.dmg
#
# The app only runs on macOS versions at least as new as the libraries it was built with,
# see the output of  vtool -show-build  and the notes in README.md.

set -e

TOPDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
BUILD=${1:-build-clean}
VERSION=${2:-$(git -C "$TOPDIR" describe --tags --match "v[0-9]*" --always 2>/dev/null || echo dev)}
XTB_ROOT=${XTB_ROOT:-$TOPDIR/../xtb-install}
QTBIN=/opt/homebrew/opt/qt@5/bin

SRC="$TOPDIR/$BUILD/bin/IQmol.app"
DIST="$TOPDIR/dist"
# Work outside the project folder: a folder that iCloud syncs (Desktop, Documents) gives the
# files attributes that code signing refuses.
WORK=$(mktemp -d /tmp/iqmolx-package.XXXXXX)
APP="$WORK/IQMolX.app"
DMG="$DIST/IQMolX-$VERSION-macOS-arm64.dmg"

[ -d "$SRC" ] || { echo "No $SRC, build IQMolX first (build-macos.sh)"; exit 1; }

rm -f "$DMG"
mkdir -p "$DIST"
ditto --norsrc --noextattr "$SRC" "$APP"

# The xtb program, so that the 'xtb program' backend works without a separate installation
if [ -x "$XTB_ROOT/bin/xtb" ]; then
   cp "$XTB_ROOT/bin/xtb" "$APP/Contents/MacOS/xtb"
else
   echo "Warning: no xtb program in $XTB_ROOT/bin, it will not be part of the app"
fi

# Qt frameworks and plugins
"$QTBIN/macdeployqt" "$APP" -libpath="$XTB_ROOT/lib" -no-strip

# Everything else that points outside the app (Homebrew, the xtb install)
python3 - "$APP" "$XTB_ROOT/lib" /opt/homebrew/opt/gcc/lib/gcc/current /opt/homebrew/lib <<'PY'
import os, re, shutil, subprocess, sys

app, extra = sys.argv[1], sys.argv[2:]
frameworks = os.path.join(app, "Contents", "Frameworks")
os.makedirs(frameworks, exist_ok=True)

def run(*cmd):
    return subprocess.run(cmd, capture_output=True, text=True)

def is_macho(path):
    if os.path.islink(path) or not os.path.isfile(path):
        return False
    with open(path, "rb") as f:
        return f.read(4) in (b"\xcf\xfa\xed\xfe", b"\xca\xfe\xba\xbe", b"\xbe\xba\xfe\xca")

def macho_files():
    for root, _, names in os.walk(app):
        for n in names:
            p = os.path.join(root, n)
            if is_macho(p):
                yield p

def deps(path):
    out = run("otool", "-L", path).stdout.splitlines()[1:]
    return [l.strip().split(" (")[0] for l in out]

def rpaths(path):
    out = run("otool", "-l", path).stdout.splitlines()
    return [out[i+2].split()[1] for i, l in enumerate(out) if "LC_RPATH" in l]

def system(dep):
    return dep.startswith("/usr/lib/") or dep.startswith("/System/")

origin = {}     # copied file -> where it came from (run paths are relative to that)

def resolve(dep, owner):
    owner = origin.get(owner, owner)
    if dep.startswith("@rpath/"):
        name = dep[len("@rpath/"):]
        for r in rpaths(owner) + extra:
            r = r.replace("@loader_path", os.path.dirname(owner)).replace(
                "@executable_path", os.path.join(app, "Contents", "MacOS"))
            if os.path.exists(os.path.join(r, name)):
                return os.path.realpath(os.path.join(r, name))
        return None
    if dep.startswith("@"):
        return None
    return os.path.realpath(dep) if os.path.exists(dep) else None

inside = lambda p: os.path.realpath(p).startswith(os.path.realpath(app) + os.sep)
done = set()
queue = list(macho_files())
while queue:
    f = queue.pop()
    if f in done:
        continue
    done.add(f)
    for dep in deps(f):
        if system(dep) or dep.startswith("@executable_path") or dep.startswith("@loader_path"):
            continue
        real = resolve(dep, f)
        if real is None:
            print("  cannot resolve", dep, "for", os.path.relpath(f, app))
            continue
        if inside(real):
            target = real
        else:
            target = os.path.join(frameworks, os.path.basename(real))
            if not os.path.exists(target):
                shutil.copy(real, target)
                origin[target] = real
                os.chmod(target, 0o755)
                run("install_name_tool", "-id", "@executable_path/../Frameworks/" + os.path.basename(real), target)
                queue.append(target)
        rel = "@executable_path/../Frameworks/" + os.path.relpath(os.path.realpath(target), os.path.realpath(frameworks))
        run("install_name_tool", "-change", dep, rel, f)

# the install names of the plugins and libraries point into the bundle, not at Homebrew
for f in macho_files():
    if f.endswith((".dylib", ".so")) and run("otool", "-D", f).stdout.splitlines()[1:]:
        run("install_name_tool", "-id", "@executable_path/../" + os.path.relpath(f, os.path.join(app, "Contents")), f)

# no run paths into the build or Homebrew trees
for f in macho_files():
    for r in rpaths(f):
        if not r.startswith("@"):
            run("install_name_tool", "-delete_rpath", r, f)
PY

# What is left that points outside the app
LEFT=$( for f in $(find "$APP" -type f \( -perm -u+x -o -name "*.dylib" -o -name "*.so" \)); do
   file -b "$f" | grep -q "Mach-O" && otool -L "$f" | tail -n +2 | awk '{print $1}' | grep -E "^/(opt|usr/local|Users)" ; done | sort -u )
if [ -n "$LEFT" ]; then echo "Still linked to libraries outside the app:"; echo "$LEFT"; exit 1; fi

xattr -cr "$APP"
# Apple Silicon only runs signed code: sign every piece ad hoc, then the app
for f in $(find "$APP" -type f); do
   file -b "$f" | grep -q "Mach-O" && codesign --force --sign - "$f" >/dev/null 2>&1 || true
done
codesign --force --deep --sign - "$APP"
codesign --verify --deep --strict "$APP" && echo "signature ok"

# Disk image with a link to /Applications
mkdir -p "$WORK/dmgroot"
cp -R "$APP" "$WORK/dmgroot/"
ln -s /Applications "$WORK/dmgroot/Applications"
hdiutil create -volname "IQMolX $VERSION" -srcfolder "$WORK/dmgroot" -ov -format UDZO "$DMG" >/dev/null
echo "The signed app is in $APP"
echo "Built: $DMG ($(du -h "$DMG" | cut -f1))"
