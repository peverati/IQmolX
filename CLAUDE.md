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
changes themselves have not run yet.  Not done: updating the helper in an existing environment
after an app update (the environment keeps the version it was made with).  TorchANI/ANI-2x was run for real (torch 2.13, TorchANI 2.9; its weights
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
Tests: `PythonEngineTest`, `PythonSetupTest`, `PythonDialogTest` (C++), `pytest` in the helper
(`IQMOLX_TEST_DOWNLOADS=1` runs the TorchANI model tests).  Never run `XtbProcessTest` in a
folder with content: it removes its working directory when xtb is not found.

## What was verified and what was not

* The macOS app builds, passes its tests and has been packaged into a signed disk image; a
  GitHub workflow (`.github/workflows/macos-release.yml`) builds it on macOS 14.
* Speed of UMA-like models on a Mac CPU/GPU, TorchANI's and AIMNet2's licenses as shipped, and
  ONNX/TorchScript export of the large models are **unverified**.
