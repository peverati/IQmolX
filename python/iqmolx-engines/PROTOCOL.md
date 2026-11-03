# IQMolX engine protocol, version 1

Copyright (C) 2026 Roberto Peverati. Developed with the assistance of Claude (Anthropic).
Licensed under the GNU General Public License v3 or later, like IQMolX.

This document describes how IQMolX talks to the Python helper (`iqmolx-engines`). The helper
runs as a child process that stays alive, so a model is loaded once and reused.

## Transport

* The helper is started as `python -m iqmolx_engines` (or the `iqmolx-engines` command).
* Messages are UTF-8 JSON objects, **one per line**, terminated by `\n`. A message never
  contains a raw newline (the helper writes compact JSON).
* IQMolX writes requests to the helper's **stdin**; the helper writes responses and events to
  its **stdout**. Nothing else is ever written to stdout (the helper redirects `print` from
  engine code to stderr). **stderr** carries free-form logs; IQMolX shows its tail on errors.
* Numbers must be finite. `NaN` and `Infinity` are not valid JSON and are never sent.
* End of input (stdin closed) makes the helper exit with status 0.

## Units and conventions

| Quantity     | Unit                                    |
|--------------|-----------------------------------------|
| coordinates  | Angstrom                                |
| energy       | eV                                      |
| forces       | eV/Angstrom (force = minus gradient)    |
| Hessian      | eV/Angstrom^2, 3N x 3N, row-major nested lists, atom-major order (x1 y1 z1 x2 ...) |
| dipole       | e*Angstrom (elementary charge times Angstrom) |
| charges      | e (partial atomic charges)              |
| charge       | integer total charge, in e              |
| multiplicity | integer 2S+1 (1 = singlet)              |

IQMolX converts to and from its own units (Hartree, Bohr) on its side; the helper stays in
the ASE-native units above so adapters need no conversion.

## Message envelope

Every request is an object with a `type` and an `id` chosen by IQMolX (an integer or string,
unique among outstanding requests):

    {"id": 7, "type": "compute", ...}

The helper answers each request with exactly one **final** message carrying the same `id`:

    {"id": 7, "ok": true,  "result": { ... }}
    {"id": 7, "ok": false, "error": {"code": "...", "message": "...", "detail": "..."}}

Before the final message the helper may send any number of **events** for that request:

    {"id": 7, "event": "progress", "fraction": 0.4, "message": "Downloading weights"}

`fraction` is optional (0 to 1; absent when unknown). Events never end a request.

The helper also sends one unsolicited message right after it starts, before any request:

    {"event": "ready", "protocol": 1, "helper": "0.1.0", "python": "3.12.7"}

Unknown fields in any message must be ignored by both sides; that is how version 1 grows.

## Versioning

`protocol` is a single integer. Version 1 only gets additive changes (new optional fields,
new request types, new error codes, new property names). A change that breaks existing
messages makes the number 2. IQMolX sends `hello` with the highest version it speaks; the
helper answers with the version it will use (the lower of the two) or the error
`unsupported_protocol` when it cannot talk any version IQMolX understands.

## Requests

### `hello`

    {"id": 1, "type": "hello", "protocol": 1, "client": "IQMolX 3.0"}

Result: `{"protocol": 1, "helper": "0.1.0", "python": "3.12.7", "platform": "darwin-arm64"}`.

### `list_engines`

    {"id": 2, "type": "list_engines"}

Result: `{"engines": [ <engine descriptor>, ... ]}`. An engine whose Python dependencies
are missing is still listed, with `"available": false` and `"unavailable_reason"`, so that
IQMolX can tell the user what to install.

#### Engine descriptor (capabilities)

The IQMolX window is built entirely from this object; IQMolX has no engine-specific code.

    {
      "name": "demo",                    // identifier used in load requests
      "title": "Demo (Morse toy potential)",
      "description": "Pairwise Morse potential for tests; not a real method.",
      "version": "0.1.0",
      "license": "GPL-3.0-or-later",     // license of the engine and its weights, free text
      "citation": "",                    // optional
      "available": true,
      "unavailable_reason": "",
      "elements": [1, 2, 3, ..., 18],    // atomic numbers; null means all elements
      "charge": true,                    // accepts a nonzero total charge
      "multiplicity": true,              // accepts multiplicity other than 1
      "properties": ["energy", "forces", "hessian", "dipole", "charges"],
      "options": [ <option>, ... ],
      "requires_download": false,        // load may fetch weights (progress events follow)
      "requires_login": false,           // gated weights (e.g. Hugging Face), see below
      "requirements": []                 // pip requirements of the engine, e.g. ["torch", "torchani"]
    }

`properties` lists what `compute` may be asked for. `energy` and `forces` are mandatory for
every engine. Without `hessian`, IQMolX builds a numerical Hessian from forces
(`Vibrations`); with it, the helper's analytic Hessian is used.

`requirements` lists what has to be installed in the Python environment for the engine to
become available (empty when the helper itself is enough). IQMolX uses it to offer the
installation of an unavailable engine: it runs the package installer of the environment on
exactly these strings and then asks for `list_engines` again. The strings are plain pip
requirement specifiers and are never executed by the helper.

An **option** describes one setting in the window:

    {"name": "well_depth", "title": "Well depth (eV)", "type": "float",
     "default": 4.0, "min": 0.1, "max": 20.0, "description": "...", "scope": "load"}

* `type`: `bool`, `int`, `float`, `string` or `choice` (with `"choices": ["a", "b"]`).
* `scope`: `load` (applied when the model loads; changing it needs a new `load`) or `compute`
  (sent with every compute request). Default `load`.
* `min`, `max`, `choices`, `description` are optional. `default` is always present.

### `load`

    {"id": 3, "type": "load", "engine": "demo", "options": {"well_depth": 4.0}}

Loads (and if needed downloads) the engine. Only one engine is active at a time; loading
replaces the previous one. Result: `{"engine": "demo", "info": {...}}` where `info` is free-form
text-friendly data (device, model name) for the status line. May emit `progress` events.
`options` holds only `load`-scope options; missing ones take their defaults.

Errors: `unknown_engine`, `engine_unavailable`, `login_required`, `bad_option`, `cancelled`,
`engine_error`.

### `compute`

    {"id": 4, "type": "compute",
     "numbers": [8, 1, 1],
     "coordinates": [[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]],
     "charge": 0, "multiplicity": 1,
     "want": ["energy", "forces"],
     "options": {}}

* `want` defaults to `["energy", "forces"]`. `options` holds `compute`-scope options.
* The helper validates the request against the loaded engine's capabilities.

Result (only the requested properties appear; arrays have the shapes below):

    {"energy": -12.3, "forces": [[fx, fy, fz], ...],         // N x 3
     "hessian": [[...], ...],                                // 3N x 3N
     "dipole": [dx, dy, dz], "charges": [q1, ..., qN]}

Errors: `not_loaded`, `bad_request`, `unsupported_element`, `unsupported_charge`,
`unsupported_multiplicity`, `unsupported_property`, `cancelled`, `engine_error`.

### `cancel`

    {"id": 5, "type": "cancel", "target": 4}

Asks the helper to abandon the request with id `target`. Semantics:

* `cancel` is answered **immediately** (`{"id": 5, "ok": true, "result": {"found": true}}`,
  `found` false when the target already finished or never existed), even while a computation
  is running; the helper reads requests on a separate thread for this reason.
* The helper sets a cancel flag. Engines that can poll it (loops, downloads, optimisation
  inside the engine) stop early. A single forward pass of a neural network cannot be
  interrupted; it runs to completion.
* In every case the target request then ends with the error `cancelled` and any result is
  discarded. IQMolX should not wait for the target to finish: after sending `cancel` it can
  mark the job as cancelled and ignore later messages with that id.
* Requests are processed **one at a time, in arrival order**. A cancelled request that is
  still queued is dropped and answered with `cancelled` without running.

### `shutdown`

    {"id": 6, "type": "shutdown"}

The helper answers `{"ok": true, "result": {}}` and exits with status 0 once the current
computation finished (or immediately when idle).

## Errors

`error.code` is a stable machine-readable string, `message` is a one-line text suitable for
the window, `detail` is optional extra text (for example a Python traceback) for a log view.

| code                     | meaning |
|--------------------------|---------|
| `bad_request`            | malformed message, missing or wrongly typed field, unknown `type`; `id` may be null when it could not be read |
| `unsupported_protocol`   | no common protocol version |
| `unknown_engine`         | no engine with that name |
| `engine_unavailable`     | Python dependencies of the engine are missing |
| `login_required`         | gated weights need a login (the user logs in from their own terminal) |
| `not_loaded`             | `compute` before a successful `load` |
| `bad_option`             | unknown option name or value of the wrong type or out of range |
| `unsupported_element`    | an atomic number the engine cannot handle (`detail` lists them) |
| `unsupported_charge`     | nonzero charge for an engine that is neutral only |
| `unsupported_multiplicity` | multiplicity other than 1 for an engine without spin |
| `unsupported_property`   | `want` contains a property not in `properties` |
| `cancelled`              | the request was cancelled |
| `engine_error`           | the engine raised an exception; the helper stays alive and usable |
| `internal`               | bug in the helper itself; the helper stays alive |

Invalid JSON on a line yields `{"id": null, "ok": false, "error": {"code": "bad_request", ...}}`
and the helper keeps reading.

### A crashed helper

When the helper process dies (killed, out of memory, a segfault in a native library), no
final message arrives. IQMolX detects the end of the stdout stream or a non-zero exit code
and reports: the exit status, the last lines of stderr, and the request that was in flight.
It does not restart the helper automatically; the user presses a button, and the model
loads again. A Python exception inside an engine is **not** a crash: it becomes an
`engine_error` response.

If the helper cannot even start (import error in the package itself), it prints a one-line
JSON `{"event": "fatal", "message": "..."}` to stdout when it still can, and exits with
status 2.

## Example session

    <- {"event":"ready","protocol":1,"helper":"0.1.0","python":"3.12.7"}
    -> {"id":1,"type":"hello","protocol":1,"client":"IQMolX 3.0"}
    <- {"id":1,"ok":true,"result":{"protocol":1,"helper":"0.1.0","python":"3.12.7","platform":"darwin-arm64"}}
    -> {"id":2,"type":"list_engines"}
    <- {"id":2,"ok":true,"result":{"engines":[{"name":"demo", ...}]}}
    -> {"id":3,"type":"load","engine":"demo","options":{}}
    <- {"id":3,"ok":true,"result":{"engine":"demo","info":{"device":"cpu"}}}
    -> {"id":4,"type":"compute","numbers":[1,1],"coordinates":[[0,0,0],[0,0,0.74]],"charge":0,"multiplicity":1,"want":["energy","forces"]}
    <- {"id":4,"ok":true,"result":{"energy":-3.99,"forces":[[0,0,0.01],[0,0,-0.01]]}}
    -> {"id":5,"type":"shutdown"}
    <- {"id":5,"ok":true,"result":{}}

## Plugins

Engines are Python classes deriving from `iqmolx_engines.base.Engine`. Built-in engines are
registered in `iqmolx_engines.registry`; third-party packages can add engines through the
entry-point group `iqmolx_engines` (each entry point names an `Engine` subclass).
