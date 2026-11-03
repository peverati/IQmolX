#!/bin/bash
# Copyright (C) 2026 Roberto Peverati.  Part of IQmol, free software under the GNU General
# Public License, version 3 or (at your option) any later version.
# Build and install xtb (Grimme group) from source on Apple Silicon macOS with Homebrew.
# IQmol uses the resulting `xtb` program (set its path in IQmol > Preferences) and, later,
# libxtb and xtb.h from the same install.
#
# Prerequisites:  brew install cmake gcc openblas
#
# OpenBLAS is used instead of Apple's Accelerate framework: with Homebrew GCC on Apple
# Silicon the Accelerate build of xtb produces NaN steps in its geometry optimizer for
# molecules of about 20 atoms and more (e.g. caffeine), while energies are still correct.
#
# Usage:  ./build-xtb-macos.sh [source-dir] [install-dir]
#   source-dir   default: <IQmol3>/../xtb
#   install-dir  default: <IQmol3>/../xtb-install
#
# The result is   <install-dir>/bin/xtb   <install-dir>/lib/libxtb.dylib   <install-dir>/include/xtb.h

set -e

TOPDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
SRC=${1:-$TOPDIR/../xtb}
PREFIX=${2:-$TOPDIR/../xtb-install}

# The xtb version IQmol was developed and tested with (6.7.1 + later fixes)
XTB_COMMIT=920a35fb512642ca7e93dac340513ea3e6267223

# Apple clang has no OpenMP, so use Homebrew GCC for both C and Fortran
GCC_BIN=$(ls /opt/homebrew/bin/gcc-[0-9]* 2>/dev/null | grep -E 'gcc-[0-9]+$' | sort -V | tail -1)
GFORTRAN_BIN=$(ls /opt/homebrew/bin/gfortran-[0-9]* 2>/dev/null | sort -V | tail -1)
if [ -z "$GCC_BIN" ] || [ -z "$GFORTRAN_BIN" ]; then
   echo "Homebrew gcc/gfortran not found.  Install with: brew install gcc" >&2
   exit 1
fi

# gfortran needs to be told where the SDK is
export SDKROOT=$(xcrun --show-sdk-path)
export CC=$GCC_BIN
export FC=$GFORTRAN_BIN

if [ ! -d "$SRC/.git" ]; then
   git clone https://github.com/grimme-lab/xtb "$SRC"
fi
cd "$SRC"
git fetch --quiet origin "$XTB_COMMIT" 2>/dev/null || git fetch --quiet origin
git checkout --quiet "$XTB_COMMIT"

# xtb downloads its dependencies (mctc-lib, tblite, dftd4, ...) while configuring.
# The policy minimum is needed because some of them predate CMake 4.
# Anaconda is ignored so that CMake does not pick up its libraries.
# A cached BLAS choice from an earlier configuration would be kept otherwise
rm -f build/CMakeCache.txt

cmake -B build -G "Unix Makefiles" \
   -DBLA_VENDOR=OpenBLAS \
   -DCMAKE_PREFIX_PATH=/opt/homebrew/opt/openblas \
   -DCMAKE_BUILD_TYPE=Release \
   -DCMAKE_INSTALL_PREFIX="$PREFIX" \
   -DWITH_TESTS=OFF \
   -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
   -DCMAKE_IGNORE_PREFIX_PATH=/opt/anaconda3

cmake --build build -j"$(sysctl -n hw.ncpu)"
cmake --install build

echo "Installed xtb to $PREFIX"
"$PREFIX/bin/xtb" --version | head -3
