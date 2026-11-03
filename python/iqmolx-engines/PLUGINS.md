# Writing an engine plugin for IQMolX

Copyright (C) 2026 Roberto Peverati. Developed with the assistance of Claude (Anthropic).
Licensed under the GNU General Public License v3 or later, like IQMolX.

An **engine** is a Python class that turns atoms into an energy and forces (and, if it can, a
Hessian). IQMolX has no code for any particular method: the *Python Engines* window is built
from what the engine says about itself, and the geometry optimizer and the frequency analysis
of IQMolX run on the forces the engine returns. So adding a method to IQMolX means writing one
class and installing it. This guide shows how, using a complete working example.

* Read [PROTOCOL.md](PROTOCOL.md) for the messages between IQMolX and the helper. You do not
  need it to write an engine, the helper handles all of it.
* The example used throughout is a real package that the tests of the helper install and run:
  [`python/examples/iqmolx-engine-lennard-jones`](../examples/iqmolx-engine-lennard-jones).
  The built-in engines are further examples: `engines/demo.py` (analytic forces and Hessian, a
  cancellable delay), `engines/ase_adapter.py` (wrapping an existing library) and
  `engines/torchani_engine.py` (lazy imports, a model that is downloaded, a Hessian from
  autograd).

## The smallest engine

```python
import numpy as np
from iqmolx_engines.base import Engine, option


class LennardJonesEngine(Engine):
    name = "lj"
    title = "Lennard-Jones (example plugin)"
    description = "Noble gas clusters with a 12-6 pair potential."
    license = "GPL-3.0-or-later"
    elements = (2, 10, 18, 36, 54)          # He, Ne, Ar, Kr, Xe
    charge = False
    multiplicity = False
    properties = ("energy", "forces")
    options = (option("epsilon_scale", "Well depth scale", "float", 1.0, min=0.1, max=10.0),)

    def load(self, options, ctx):
        self.scale = options["epsilon_scale"]
        return {"model": "Lennard-Jones"}

    def compute(self, numbers, coordinates, charge, multiplicity, want, options, ctx):
        ...                                  # energy in eV, forces in eV/Angstrom
        return {"energy": energy, "forces": forces}
```

The class attributes say what the engine is; `load` prepares it; `compute` does the work. The
complete file is
[`iqmolx_engine_lj/__init__.py`](../examples/iqmolx-engine-lennard-jones/src/iqmolx_engine_lj/__init__.py).

## What the attributes mean

Everything here is sent to IQMolX as the engine's **descriptor** (`list_engines`), and the window
is built from it. `iqmolx-engines --describe` prints the descriptors of all engines as JSON.

| Attribute | Meaning |
|---|---|
| `name` | Identifier used in `load` requests and in the preferences. Short, lower case, **unique**; required. |
| `title` | What the window shows, and what the result is labelled with. |
| `description` | A sentence or two for the window. Say what the method is and what it cannot do. |
| `version` | Version of the engine (default `"0.1.0"`). |
| `license` | Free text: the licence of the engine **and of its weights**. Shown to the user. |
| `citation` | What to cite when the method is used. Shown in the window. |
| `elements` | Tuple of atomic numbers the engine handles, or `None` for all. |
| `charge` | `True` if a nonzero total charge is allowed. |
| `multiplicity` | `True` if a spin multiplicity other than 1 is allowed. |
| `properties` | What `compute` can return: any of `"energy"`, `"forces"`, `"hessian"`, `"dipole"`, `"charges"`. `energy` and `forces` are mandatory. |
| `options` | The settings, see below. |
| `requires_download` | `True` if loading may fetch weights: the window tells the user. |
| `requires_login` | `True` for gated weights (see *Weights and licences*). |
| `requirements` | pip requirement strings the engine needs, for example `("torch", "torchani")`. When the engine is unavailable IQMolX offers an *Install* button that installs exactly these. |

The helper **checks the request against these attributes before your code runs**: an element,
a charge, a spin or a property you did not declare is refused with `unsupported_element`,
`unsupported_charge`, `unsupported_multiplicity` or `unsupported_property`, so `compute` can
assume a valid request. (This is `Engine.validate_system`, you can call it in tests.)

### Charge and spin

`charge` and `multiplicity` reach `compute` as plain integers; passing them on is the engine's job.
If you wrap a library whose calculator takes them as arguments, set them there on every call, and
reset anything the library caches (the ASE adapter, `engines/ase_adapter.py`, shows both). Engines
that cannot handle them declare `charge = False` and `multiplicity = False`, and the helper refuses
such requests for you, so a charged molecule is never silently calculated as a neutral one.

### Properties and what IQMolX does with them

* `forces` drive the geometry optimization (IQMolX's own L-BFGS) and, if there is no Hessian,
  the vibrational frequencies by finite differences of forces. That works, but the step of
  0.01 Bohr costs 6N force evaluations and is not as accurate as an analytic Hessian: for ANI-2x
  it changes the largest Hessian elements by a few tenths of an eV/Angstrom squared.
* `hessian` is used directly when you offer it. If the library can give you one (autograd,
  analytic), do: `torchani_engine.py` does it with `torch.autograd.functional.hessian`.
* `dipole` and `charges` appear in the result summary.

### Options

`option(name, title, type, default, scope="load", **extra)` describes one setting. The window
makes a widget for each: `bool` a check box, `int` and `float` a spin box (with `min` and `max`),
`string` a text field, `choice` a menu (with `choices=[...]`). `description=` becomes the tooltip.

* `scope="load"`: applied when the model loads; changing it loads the model again.
* `scope="compute"`: sent with every calculation (`options` argument of `compute`).
* **Names must be unique within an engine**, whatever the scope: the window stores the values
  by name.
* Values are checked before your code runs (type, range, choice); a bad value is the error
  `bad_option`. Defaults are filled in, so `options` in `load` and `compute` is always complete.
* An option that has no sensible default (a model name, a class to use) is declared with
  `required=True` and a `placeholder="an example value"`: the window then shows the hint in the
  empty field and keeps *Run* off until it is filled in, and the helper refuses an empty value for you
  (`bad_option`, "'Title' is required"). The ASE adapter does this for its calculator class.

## `load` and `compute`

```python
def load(self, options, ctx):            # options: the "load" options, complete
    ...
    return {"device": "cpu", "model": "..."}     # shown in the progress text; may be None

def compute(self, numbers, coordinates, charge, multiplicity, want, options, ctx):
    ...
    return {"energy": e, "forces": f}
```

* The helper makes **one instance** per `load` request (`cls()`, so `__init__` takes no
  arguments) and keeps it until the next `load`. Keep state in attributes set in `load`.
* `numbers` is an integer numpy array (N,), `coordinates` a float array (N, 3) in **Angstrom**,
  `charge` and `multiplicity` are ints, `want` the list of property names asked for.
* Return a dict with **at least every property in `want`**; extra keys are ignored. Anything
  missing is the error `engine_error`.

| Quantity | Unit and shape |
|---|---|
| `energy` | eV, a number |
| `forces` | eV/Angstrom, (N, 3); force = minus gradient |
| `hessian` | eV/Angstrom^2, (3N, 3N), atom-major (x1 y1 z1 x2 ...) |
| `dipole` | e*Angstrom, (3,) |
| `charges` | e, (N,) |

The helper checks the shapes and that every number is finite; a wrong shape or a `nan` is
reported as `engine_error` with a message that names the property. Convert from your library's
units yourself (Hartree to eV is 27.211386245988).

### Progress, cancel and errors: `ctx`

`ctx` is an `iqmolx_engines.base.Context`:

* `ctx.progress(fraction, message)` shows progress in the window (`fraction` between 0 and 1,
  or `None`). Use it for anything that takes more than a moment: downloading weights, loading
  a large model.
* `ctx.check_cancel()` raises the cancel error if the user pressed Cancel; `ctx.cancelled()`
  just tells you. Call it in loops and between the steps of a download. **A single call into a
  library (a neural network forward pass) cannot be interrupted**: it runs to its end, and the
  helper then discards the result and answers "cancelled" (see *Cancel* in PROTOCOL.md).

Errors:

* To refuse something, raise `iqmolx_engines.base.EngineError(code, message, detail="")`. Use a
  code of the protocol: `bad_option`, `login_required`, `engine_unavailable`, `engine_error`, and so
  on (the table is in PROTOCOL.md). `message` is one line for the window, `detail` may hold more
  (a traceback, a hint).
* Any other exception that escapes your code becomes `engine_error` with the traceback in
  `detail`. **The helper stays alive**: an engine that fails does not take the program down.
* A crash of native code (a segfault, running out of memory) does kill the helper process;
  IQMolX then reports the exit status and the end of stderr.
* **Never write to stdout.** It carries the protocol. The helper points both Python's `print` and
  the file descriptor that native libraries write to (a C or Fortran library printing its own
  progress) at stderr for you, and stderr is shown to the user when something goes wrong, but use
  `ctx.progress` for anything the user should read. A library that prints a lot can usually be
  silenced with an option (TBLite has `verbosity`).

### An engine whose packages may be missing

`unavailable_reason()` is a class method that returns `""` when the engine can run and a sentence
otherwise. It is called every time the engines are listed, so **do not import heavy packages
there**: look for them with `importlib.util.find_spec`, import them in `load`.

```python
@classmethod
def unavailable_reason(cls):
    missing = [m for m in ("torch", "torchani") if importlib.util.find_spec(m) is None]
    return "%s not installed" % " and ".join(missing) if missing else ""
```

An engine that is not available is still listed. The window shows the reason, and, with
`requirements = ("torch", "torchani")`, an *Install torch, torchani* button that runs
`uv pip install torch torchani` in the environment (offered for a virtual environment, such as the
one the *Set Up Environment* button makes, when `uv` is available).

## Packaging a plugin

A plugin is an ordinary Python package that registers its class in the entry-point group
`iqmolx_engines`. The whole `pyproject.toml` of the example:

```toml
[build-system]
requires = ["setuptools>=61"]
build-backend = "setuptools.build_meta"

[project]
name = "iqmolx-engine-lennard-jones"
version = "0.1.0"
requires-python = ">=3.9"
license = { text = "GPL-3.0-or-later" }
dependencies = []

[project.entry-points."iqmolx_engines"]
lj = "iqmolx_engine_lj:LennardJonesEngine"

[tool.setuptools.packages.find]
where = ["src"]
```

* The value on the right is `module:Class`; the name on the left is free. The engine's identity
  is the `name` attribute of the class.
* **Do not list `iqmolx-engines` as a dependency.** The helper is already in the environment
  that runs the plugin (IQMolX installs it there from a wheel in the application, it is not on
  PyPI), so a requirement on it would make the installer look for it and fail. List only what
  your engine adds. `numpy` comes with the helper.
* Install the plugin **into the environment IQMolX uses** (the *Python* field of the window; for
  the environment IQMolX made, `~/Library/Application Support/IQmol/python`):

  ```bash
  uv pip install --python "$HOME/Library/Application Support/IQmol/python/bin/python" ./my-plugin
  ```

  or, while developing, `-e ./my-plugin` so that edits are picked up.
* Press **Find Engines** in the window: a running helper notices a plugin that was installed
  after it started. **Code of a plugin that was already imported is not reloaded**: after changing
  it, restart IQMolX (the helper stops with it).

### What the helper does with plugins

At every listing it loads all entry points of the group. A plugin that cannot be used is skipped
with a line on stderr that names it, and the other engines still work:

| Problem | Message (stderr) |
|---|---|
| import fails | `cannot load plugin <name>: <error>` |
| not a subclass of `Engine` | `plugin <name> is not an Engine subclass, ignored` |
| empty `name` attribute | `plugin <name> has no engine name, ignored` |
| same `name` as a built-in engine | `the plugin engine '<name>' is hidden by the built-in engine of that name` |

Built-in names (`demo`, `ase`, `torchani`) win. Pick a different name.

## Testing a plugin

1. **Look at the descriptor**: `iqmolx-engines --describe` (or
   `python -m iqmolx_engines --describe`, with the interpreter of the environment). Your engine
   must be listed with `"available": true` and the attributes you meant.
2. **Test the engine in a normal pytest**: build it, call `load` and `compute` with a
   `Context()`. The sample's tests
   ([`test_plugin_example.py`](tests/test_plugin_example.py)) show the usual ones: a case with a
   known answer, **forces against a finite difference of the energy**, the Hessian the same way
   if you offer one, and `validate_system` for what you refuse. `Context(cancelled=lambda: True)`
   tests cancellation.

   ```python
   eng = make()
   f = eng.compute(z, xyz, 0, 1, ["energy", "forces"], {}, Context())["forces"]
   ```
3. **Run it through the real helper**, as IQMolX does: one JSON object per line on stdin, one
   per line back (PROTOCOL.md has a whole session). `test_plugin_example.py` also does this for
   an installed plugin.

   ```bash
   printf '%s\n' '{"id":1,"type":"load","engine":"lj"}' \
     '{"id":2,"type":"compute","numbers":[18,18],"coordinates":[[0,0,0],[0,0,3.8]]}' \
     | python -m iqmolx_engines
   ```
4. **Try it in IQMolX**: *Calculation > Python Engines*, *Find Engines*, choose the engine,
   optimize a small molecule. A good force implementation converges in a few cycles; if the
   optimizer wanders, the forces are not the gradient of the energy (test 2 shows it).

## Weights and licences

* **Never put weights or big libraries in the plugin or in IQMolX.** Download them when the
  engine loads, tell the user (`requires_download = True`, and `ctx.progress` while downloading,
  with `ctx.check_cancel()` between chunks) and let the library cache them. TorchANI is the
  model for this: its package is 0.5 MB and it fetches the 55 MB parameter file from Hugging Face
  the first time.
* **Say what the licence is**, of the code and of the weights, in `license`, and the paper in
  `citation`. Check the weights, not only the code: they often have their own licence.
* **Gated weights** (a login or an agreement is needed): set `requires_login = True` and, when
  the model cannot be fetched, raise `EngineError("login_required", "...")` that tells the user to log
  in **in their own terminal** (for example `huggingface-cli login`). Never ask for a password or a
  token in IQMolX, never write one to a file, a preference or the output.
* Keep the plugin's own licence compatible with how you distribute it. IQMolX is GPL v3 or later;
  a plugin is a separate package and may have another licence, but it is installed next to
  GPL code and should say so clearly.

## A checklist

- [ ] `name` is unique and not one of the built-in names.
- [ ] `elements`, `charge` and `multiplicity` say exactly what the method supports.
- [ ] `properties` lists only what `compute` really returns; `forces` are minus the gradient.
- [ ] Units: eV, eV/Angstrom, Angstrom in; shapes as in the table.
- [ ] Forces agree with finite differences of the energy (and the Hessian with those of the forces).
- [ ] `unavailable_reason()` imports nothing heavy; `requirements` lists what to install.
- [ ] Long steps call `ctx.progress` and `ctx.check_cancel`.
- [ ] Nothing is written to stdout; errors use `EngineError` with a code.
- [ ] `license` and `citation` are filled in; weights are downloaded, not bundled.
- [ ] `iqmolx-engines --describe` shows the engine, and it runs through the helper and in the window.
