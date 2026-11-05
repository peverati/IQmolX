# IQMolX

IQMolX is a fork of [IQmol 3](https://github.com/nutjunkie/IQmol3), the molecular
builder and visualization package written by Andrew Gilbert.  It is called IQMolX to
tell it apart from the original program.

## How IQMolX differs from IQmol

Everything IQmol does is still there (building molecules, setting up and submitting
Q-Chem jobs, surfaces, orbitals, vibrations, ...).  IQMolX adds semiempirical
calculations with [xtb](https://github.com/grimme-lab/xtb) (Grimme group), run
from inside the program:

* **xtb window** (Calculation → xtb Setup, or the xtb button on the tool bar) to set
  up and run GFN2-xTB and GFN1-xTB calculations: single point energies, geometry
  optimizations, vibrational frequencies and optimization plus frequencies.  Charge and
  spin multiplicity, the ALPB implicit solvents and the electronic temperature can be
  chosen, and the live xtb output and a summary of the result are shown in the window.
* **Two backends**: the built-in library (libxtb, with its own geometry optimizer; the
  default) or the external xtb program.  The program adds the electronic temperature, the
  molecular orbitals, extra command line options and a working directory, and *Save All
  Files* copies every file xtb wrote.
* **Results in the usual tools**: the optimized geometry (undoable), the vibrational
  modes and IR intensities (in the Frequencies viewer) and the molecular orbitals (in
  the orbital viewer, with the program backend) are added to the molecule.
* **g-xTB (development version)** as a third backend of the xtb window, for the new
  [g-xTB](https://github.com/grimme-lab/g-xtb) method of the Grimme group: a Download button
  fetches the g-xTB build of xtb (about 41 MB, checked against its SHA-256 sum) the first time.
  It is a preliminary version and is kept apart from the normal xtb backends.
* **GFN-FF as a force field**, next to UFF and the others in Build → Select Force Field.
* **Python engines window** (Calculation → Python Engines, Ctrl+Shift+E) to calculate with
  machine learning potentials and other Python methods through one small helper
  (`python/iqmolx-engines`): the TorchANI models (ANI-1x, ANI-1ccx, ANI-2x and others), the
  Egret-1 models, Hartree-Fock and DFT with PySCF (including wB97M-V and any functional you
  name), and any ASE calculator.  Energies, geometry optimization and frequencies work like
  in the xtb window.  Third parties can add engines as plugins, see
  `python/iqmolx-engines/PLUGINS.md`.
* **Apple Silicon (arm64) macOS builds** with current Homebrew and Xcode:
  `build-macos.sh` and, for xtb, `build-xtb-macos.sh`.
* Smaller changes: visible check boxes in the Model View on macOS dark mode, clicking
  one of several orbital surfaces shows it alone, a readable (white) Q-Chem input
  preview in dark mode, and the def2-QZVP, def2-QZVPP, def2-QZVPD and def2-QZVPPD basis
  sets in the Q-Chem basis list.

<p align="center"><img src="resources/help/images/XtbProgram.png" width="420" alt="The xtb window with the xtb program as backend"></p>

The xtb help page is part of the built-in help (Help → Show Help).  The changes are on
top of IQmol 3.2.x; internally the program, its settings and its source folders keep the
name IQmol.

## Download (macOS, Apple Silicon)

Ready-made apps are on the [Releases page](https://github.com/peverati/IQmolX/releases).
They need a Mac with Apple Silicon (M1 or later) and macOS 15 or newer.  Intel Macs, Windows
and Linux have no download yet: build from source (below).

1. Download the `IQMolX-...-macOS-arm64.dmg` file of the latest release, open it and drag
   IQMolX to the Applications folder.  Everything it needs is inside, including xtb.
2. The first time, macOS will say it cannot check the app, because it is not notarized by
   Apple.  Right-click (or Control-click) IQMolX in Applications, choose **Open**, and
   confirm with **Open**.  If that does not work, run this once in a terminal:
   `xattr -dr com.apple.quarantine /Applications/IQMolX.app`
3. Start with Calculation → xtb Setup, or the xtb button on the tool bar.

### Python engines (needs the internet the first time)

xtb works right away.  The Python engines need a Python environment, which is not part of the
download (PyTorch alone is hundreds of megabytes).  Open Calculation → Python Engines and press
**Set Up Environment**: the program `uv` that is inside the app downloads its own Python and
makes an environment in `~/Library/Application Support/IQmol/python`; no Python has to be
installed.  An engine that needs more packages (TorchANI, Egret, PySCF) has an **Install**
button, and the model weights of an engine are downloaded the first time it is used.  You can
also point the window at your own Python or conda environment.  Delete the folder above to
remove everything.  Step by step, with the licenses of the engines: Help → Show Help →
Python engines.

Problems and questions are welcome as [issues](https://github.com/peverati/IQmolX/issues).

## How IQMolX was developed

The changes to IQmol were written by Roberto Peverati **with the help of Claude**, an AI
assistant made by Anthropic (Claude Code), which wrote most of the code, tests and
documentation under his direction.  He defined the goals, ran and checked the program, and
decided what to keep.  The work is tested (see `tests`), but it has not had the review of a
long-established project, so please report problems.  The commit messages of the AI-assisted
changes end with a `Co-Authored-By: Claude` line.

## Building

The source relies on submodules, so use the recursive flag to check out the code:

```
git clone --recursive https://github.com/peverati/IQMolX.git
```

To compile, make sure that your Qt installation can be found by cmake.  This means that
the CMAKE\_PREFIX\_PATH environment variable should include the directory containing the
Qt5Config.cmake file

```
export CMAKE_PREFIX_PATH=/directory/containing_Qt5Config.cmake
./configure
cd build
make
```

### macOS on Apple Silicon

`build-macos.sh` holds a recipe that works with Homebrew (qt@5, boost, openssl@3, libomp,
gcc; the header of the script lists the prerequisites).  To get the built-in xtb library,
first build xtb with `build-xtb-macos.sh` and give its install folder to the IQMolX build:

```
./build-xtb-macos.sh
XTB_ROOT=<xtb install dir> ./build-macos.sh
```

The external xtb program can also be any xtb 6.7 installation; its location is set in the
xtb window.  The tests of the xtb interface are in `tests`.

The Python helper is tested separately: `pip install -e "python/iqmolx-engines[test]" ase`
and `pytest python/iqmolx-engines` (set `IQMOLX_TEST_DOWNLOADS=1` to also run the
model tests, which download weights).  `package-macos.sh` bundles the helper as a wheel and
`uv`, see the header of the script.

## Credits and license

IQmol was written by Andrew Gilbert, see <http://iqmol.org> for the original program and
pre-compiled binaries; a user guide is in the `doc` directory.  The xtb additions
are by Roberto Peverati.  IQMolX uses
[xtb](https://github.com/grimme-lab/xtb) (LGPL-3.0) and, like IQmol, is distributed under
the GNU General Public License, see `LICENSE`.
