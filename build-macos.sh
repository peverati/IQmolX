#!/bin/bash
# Copyright (C) 2026 Roberto Peverati.  Part of IQmol, free software under the GNU General
# Public License, version 3 or (at your option) any later version.  Developed with the
# assistance of Claude (Anthropic).
# Build IQmol on Apple Silicon macOS with Homebrew (qt@5, boost, openssl@3, libomp, gcc).
#
# Prerequisites:
#   brew install cmake qt@5 boost openssl@3 libomp gcc
#   CMake 3.x in a virtualenv (the openbabel submodule is not compatible with CMake 4):
#     python3 -m venv <dir>/.venv-cmake3 && <dir>/.venv-cmake3/bin/pip install "cmake==3.31.*"
#   and point CMAKE3 at its cmake, e.g. CMAKE3=$HOME/.../.venv-cmake3/bin/cmake
#
# Usage:  ./build-macos.sh [build-dir]      (default: build-clean)
#         XTB_ROOT=<xtb install dir> ./build-macos.sh   also enables the libxtb tests

set -e

TOPDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
BUILD=${1:-build-clean}
CMAKE3=${CMAKE3:-$TOPDIR/../.venv-cmake3/bin/cmake}
BREW=/opt/homebrew/opt

# gfortran needs to be told where the SDK is
export SDKROOT=$(xcrun --show-sdk-path)

mkdir -p "$TOPDIR/$BUILD"
cd "$TOPDIR/$BUILD"
rm -f CMakeCache.txt

# Note: do not add the bare /opt/homebrew prefix, its include/ holds Qt6 headers.
"$CMAKE3" -G "Unix Makefiles" \
   -DBUILD_SHARED=OFF -DQARCHIVE=OFF -DCMAKE_BUILD_TYPE=Release \
   -DCMAKE_PREFIX_PATH="$BREW/qt@5;$BREW/libomp;$BREW/openssl@3;$BREW/boost" \
   -DCMAKE_IGNORE_PREFIX_PATH=/opt/anaconda3 \
   -DOPENSSL_ROOT_DIR=$BREW/openssl@3 \
   -DPTHREAD_LIBRARY=pthread \
   ${XTB_ROOT:+-DXTB_ROOT="$XTB_ROOT"} \
   "$TOPDIR"

make -j"$(sysctl -n hw.ncpu)" IQmol
echo "Built: $TOPDIR/$BUILD/bin/IQmol.app"
