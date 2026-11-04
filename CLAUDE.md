# IQMolX: notes for Claude Code

IQMolX is a fork of IQmol 3 (https://github.com/nutjunkie/IQmol3, by Andrew Gilbert) with xtb
support, owned by Roberto Peverati (GitHub: peverati, fork: https://github.com/peverati/IQmolX).
It is GPL v3 or later.  Qt 5 / C++17 desktop application (macOS Apple Silicon is the main
target), built with CMake.  Read README.md first.

## Rules that always apply

* **The macOS version is the main product and must never break.**  Every shared build-file
  change has to be checked with a macOS rebuild and the xtb tests.  Windows- or Linux-only
  code is guarded (`if(WIN32)`, `#ifdef Q_OS_WIN`).
* **Be upfront about AI use.**  New source files carry this header (copy it from
  `src/Util/XtbProcess.C`): `Copyright (C) 2026 Roberto Peverati`, `Developed with the
  assistance of Claude (Anthropic).` and the full GPL notice.  Commit messages end with
  `Co-Authored-By: Claude ...` (the exact line comes from the session's attribution reminder).
* Never put credentials (GitHub or Hugging Face tokens) in files, chat or preferences.  The
  user enters them in their own terminal.
* Ask before anything outward-facing: pushing, tagging, publishing releases, deleting.
* Match the surrounding code's style, comment density and naming.
* The user does not have a Windows PC.  Windows can only be checked on GitHub's runner.

## Branches and folders

* `xtb-library` (default branch on the fork): the main line.  Builds the macOS release.
* `windows-port`: experimental Windows build (MSYS2); worked on in a different chat.  Do not
  touch or merge it from here.
* `ml-engines` (this worktree, folder `IQmolX-uma` next to `IQmol3`): the Python engines work.
  Branch from `xtb-library`, merge back only after a macOS check.
* Submodules in `modules/` are not initialised in this worktree.  Run
  `git submodule update --init --recursive` before the first build.
* The xtb installation used by the builds is `../xtb-install` (sibling folder; made by
  `build-xtb-macos.sh`).  Build folder convention: `build-clean` (ignored by git).

## Build and test (macOS)

    XTB_ROOT=../xtb-install ./build-macos.sh          # needs CMAKE3 (cmake 3.31 venv), see script
    cd build-clean && make XtbProcessTest XtbBackendTest XtbDialogTest ...
    ./package-macos.sh                                # portable app + dmg (ad hoc signed)

A first build takes about an hour.  The Python helper work below does not need a build.

## The existing xtb architecture (do not change it for this work)

* `src/Util/Calculator.h`: `Calc::Input` (atomic numbers, coordinates in Angstrom, charge,
  multiplicity, solvent, electronic temperature, orbitals flag), `Calc::Result`, and the
  abstract `Calc::Calculator` with `energy()`, `optimize()`, `frequencies()`, progress callback,
  `cancel()`.
* `XtbProcess` (runs the xtb program), `XtbLibrary` (libxtb; own optimizer), shared parsing in
  `XtbOutput`, `Optimizer.h` (L-BFGS), `Vibrations.{h,C}` (mass weighting, projection, Jacobi
  eigen solver, frequencies and modes from a Hessian).
* `src/Util/XtbResultDialog.{h,C}`: the xtb window (singleton).  `Layer::Molecule::runXtb()` in
  `src/Layer/MoleculeLayer.C` runs a job and applies results (undoable geometry command, the
  Frequencies layer, orbitals).  `Preferences.{h,C}` holds settings.
* Tests are in `tests/` (standalone executables, offscreen Qt for dialogs).
* Lessons learned: range-for over a temporary `QStringList() << ...` is a dangling-reference
  bug; a Qt window must be hidden before its modality is changed; the project folder is under
  iCloud-synced Desktop, so packaging is done in `/tmp`.

## The task of this branch: a flexible Python engine helper

Goal: let IQMolX use Python-based calculators (machine learning potentials and other methods)
through one small helper, separate from xtb and from the xtb window.

Decisions already made with the user:

* xtb stays native (library/program) and keeps its own window.  The Python engines get their
  own helper, `Calc::PythonEngine` class and window (menu: Calculation; working name "Python
  engines" or "Machine learning potentials").
* The helper lives in this repository under `python/iqmolx-engines/` (installable package,
  command-line entry point).  It talks to IQMolX with **one JSON message per line over
  stdin/stdout**.  The process stays alive so a model loads once.
* Messages: list engines and their **capabilities** (supported elements, takes charge and spin,
  forces, Hessian, dipole, options with types and defaults); load a model; compute (numbers,
  coordinates, charge, spin, what to return); cancel; shut down.  The IQMolX window is built
  from the capabilities, so IQMolX contains no engine-specific code.
* Reuse the existing optimizer and `Vibrations` (numerical Hessian from forces); do not copy them.
* Plugin design: a generic ASE-calculator adapter first, then specific engines.
* **Order of work (user's choice):** 1) protocol v1 written down, 2) helper package with a
  dependency-free **demo engine** (toy potential, numpy only, for tests and CI) and the ASE
  adapter, tested from the command line, 3) IQMolX side: `Calc::PythonEngine`, capability-driven
  window, tests against the demo engine, 4) the window asks for the Python executable or
  environment (like the xtb path; automatic environment setup comes later), 5) first real engine:
  **TorchANI / ANI-2x** (H,C,N,O,F,S,Cl, neutral only; verify its license and that the
  parameters ship with the package), 6) later engines: UMA (Meta fairchem, gated weights under the
  FAIR Chemistry License, needs a Hugging Face login), AIMNet2, MACE-OFF, Egret-1, AXEL.
* Do not release the window to students before step 5 **and the user guide below**; until then it
  stays on this branch.
* Weights and PyTorch are never bundled in the app; the setup step downloads them.
* **Two guides are part of the plan** (user's request), written *as* the helper and the plugin
  discovery are built, not afterwards: a step is not finished until its part of the guides is.
  Every new engine (step 6) adds its section to both.
  1. **User guide, setting up Python** (for students): an in-app Help page like `xtb.html`
     (`resources/help/python-engines.html`, registered in `resources/IQmol.qrc`, linked from
     `resources/help/index.html` and `toolbar.html`, with screenshots of the window).  It covers
     the Python Engines window, Set Up Environment (what uv downloads, where the environment lives,
     how to remove it), the Install buttons, using your own Python or conda environment, the first
     use of an engine (weights downloaded: sizes, internet needed, proxies), the ASE adapter with a
     worked example, licences and citations of each engine, and troubleshooting (uv missing, helper
     stopped, engine not available, out of memory).
  2. **Developer guide, writing engine plugins** (`python/iqmolx-engines/PLUGINS.md`, next to
     `PROTOCOL.md`): the `Engine` base class and what each attribute does (elements, charge,
     multiplicity, properties, options, requirements, requires_download and requires_login), the
     `load` and `compute` contract and units, the error codes, progress and cancel through
     `Context`, plugin discovery through the `iqmolx_engines` entry-point group with a complete
     `pyproject.toml`, testing a plugin (`--describe`, a protocol session, the helper tests),
     licences and weights (declared, never bundled), and a checklist.  The examples are real code:
     a small sample plugin package under `python/examples/` that the helper tests install and run,
     so the guide cannot drift from the code.

Open questions to settle early: exact protocol fields and versioning, cancel semantics while a
model is computing, how a crashed helper is reported, where the Python environment lives (the
app-support folder), and whether `uv` can create it without a system Python.

**Status (steps 1 to 5 are done on `ml-engines`).**  Answers: the protocol is in
`python/iqmolx-engines/PROTOCOL.md` (v1; cancel is answered at once and a forward pass runs to
its end; a lost helper is reported with its exit status and stderr tail); the environment is
`~/Library/Application Support/IQmol/python`, made by `Calc::PythonSetup` with `uv venv --python
3.12 --python-preference only-managed`, so no system Python is needed (verified with a real uv).
`uv` itself is found in `$IQMOLX_UV`, the app's `Resources/bin`, the PATH or the usual install
folders.  `package-macos.sh` now bundles uv (from `UV_BIN`, with its license texts) and the helper
as a wheel in `Resources/python` (installing from a folder would write into the signed app); the
workflow installs a pinned uv from PyPI for that.  A packaged app was checked in a real run: with
no environment variables it built its environment from the bundled uv and wheel.  The CI workflow
changes themselves have not run yet.  Updating the helper after an app update: the window shows "Update Helper"
(`PythonSetup::updateHelper`, uv `--reinstall-package`) when the bundled helper's version is higher
than the running one (`hello` -> `helper`); only for the IQMolX-made environment.  **So bump
`__version__` in `iqmolx_engines/__init__.py` and `version` in `pyproject.toml` together whenever
the helper changes in a way users should get** (new engine, fix); without the bump the button never appears.  TorchANI/ANI-2x was run for real (torch 2.13, TorchANI 2.9; its weights
are downloaded from Hugging Face on first use, MIT), not frequencies in the app window.  The
ASE adapter and the "Install ..." button were checked with real uv outside the window's click.
**Guides.**  The developer guide `python/iqmolx-engines/PLUGINS.md` is written, with the sample
plugin `python/examples/iqmolx-engine-lennard-jones` that the helper tests install with pip and
run (this is also where the entry-point discovery was first tried with a real plugin: it works,
and a running helper notices a plugin installed after it started).  The user guide is the Help page
`resources/help/python-engines.html` (screenshots of the real window; checked by rendering it with
QTextBrowser, not seen in the app's Help window).  Each later engine adds its section to both
guides.  The release gate in the plan ("step 5 and the user guide") is now met for the guide; the
open items before a release are the ones listed above (CI run, helper update after an app update,
ANI-2x frequencies tried in the window).
**Engines added after step 5 (user's order: TorchANI models, Egret, PySCF).**  All are in
`python/iqmolx-engines/src/iqmolx_engines/engines/`.  *TorchANI* (`torchani_engine.py`): one engine
per model on a shared base (ANI-2x keeps the name `torchani`; ani1x, ani1ccx, ani1xnr, ani2xr, ani2dr,
anir2s with a solvent option, animbis with charges); all weights are MIT on Hugging Face except
`snn-ani2xr`, whose repository states no license, so it is left out.  *Egret* (`egret_engine.py`, five
engines egret1, egret1e, egret1t, egret1m, egret1s): licenses read from the primary sources (repository
LICENSE, README, preprint: MIT; mace-torch MIT), weights are single files downloaded by
`download.py` from a pinned commit and checked against the git blob hash; elements from the preprint.
*PySCF* (`pyscf_engine.py`): HF and eight DFT functionals; **PySCF 2.14's analytic Hessian is wrong for
meta-GGAs (TPSS, M06-2X, r2SCAN), so those use finite differences of the analytic forces**.  What was
run for real (user approved the downloads): all eight TorchANI engines (ANI-1xnr's energies are on a
valence scale, water -468 eV, documented; ANI-r2s only vacuum and water, chloroform and acetonitrile not
run), all five Egret models (weights downloaded and hash-checked; **they have no Hessian: the compiled
TorchScript files cannot run that routine**, so `properties` is energy and forces and IQMolX makes
frequencies from forces; element lists confirmed from the model files), PySCF, and the download helper
against a local server.  Egret and ANI need `mace-torch` / `torch`+`torchani` in the environment (Install
buttons).  The user's environment keeps
the helper it was made with: reinstall it to get new engines (see the earlier note).
Tests: `PythonEngineTest`, `PythonSetupTest`, `PythonDialogTest` (C++), `pytest` in the helper
(`IQMOLX_TEST_DOWNLOADS=1` runs the TorchANI model tests).  Never run `XtbProcessTest` in a
folder with content: it removes its working directory when xtb is not found.

## g-xTB backend (development version)

A third backend of the xtb window (`Preferences::XtbBackend() == "gxtb"`), added because g-xTB
(github.com/grimme-lab/g-xtb) is a modified xtb 6.7.1 build run with `--gxtb`, preliminary and
changing; it never replaces the normal backends.  `Calc::XtbProcess::GXTB` (command line `--gxtb`;
solvent as `--gbe`/`--cosmo`, ALPB does not work with it; no `--etemp`; **no `--uhf` for singlets**,
which would force an unrestricted wavefunction; **`OMP_NUM_THREADS=1` for Hessians**, the macOS build's
threaded numerical Hessian is wrong).  `Calc::GxtbSetup` (Util/GxtbSetup.{h,C}) is the Download
button: curl + SHA-256 check + tar into `AppDataLocation/g-xtb`, pinned to a repository commit and the
v2.0.1 build (14 May 2026, macOS arm64 and x86_64 only); to move to a newer g-xTB change the commit,
the file names and the sums in `GxtbSetup.C` together, rerun `GxtbTest` with a real program and
`XtbOrbitalsTest <g-xtb xtb> gxtb`.  Not bundled in the dmg: the license of the parameters is not stated.
The Molden file is written by tblite: its Cartesian d coefficients (xy, xz, yz) are sqrt(3) larger than
IQmol's, handled by `parseMolden(..., tbliteConvention)`.  Tests: `GxtbTest` (fake archive, command line,
threads; with a g-xtb xtb as argument also a real run; `GXTB_REAL_DOWNLOAD=1` downloads from GitHub),
`XtbOrbitalsTest <xtb> gxtb`, `XtbDialogTest`.  Windows/Linux: no download yet (macOS only), a
user-chosen g-xtb program works anywhere.

## What was verified and what was not

* The macOS app builds, passes its tests and has been packaged into a signed disk image; a
  GitHub workflow (`.github/workflows/macos-release.yml`) builds it on macOS 14.
* Speed of UMA-like models on a Mac CPU/GPU, TorchANI's and AIMNet2's licenses as shipped, and
  ONNX/TorchScript export of the large models are **unverified**.
