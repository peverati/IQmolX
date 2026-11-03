# iqmolx-engines

Helper process that lets IQMolX use Python-based calculators (machine learning potentials
and other methods). IQMolX talks to it with one JSON message per line over stdin/stdout;
see [PROTOCOL.md](PROTOCOL.md).

    pip install -e '.[test]'
    iqmolx-engines --describe        # list engines and their capabilities
    pytest                           # tests (the ASE adapter tests skip without ASE)

Engines: `demo` (Morse toy potential, numpy only, for tests), `torchani` (ANI-2x, needs torch and
torchani; the parameters are downloaded on first use) and `ase` (generic adapter
for any ASE calculator). Weights and PyTorch are never bundled with IQMolX.

To add a method, write an engine plugin: see [PLUGINS.md](PLUGINS.md) (the guide for developers,
with a complete sample plugin in [`../examples`](../examples)). The messages between IQMolX and
the helper are in [PROTOCOL.md](PROTOCOL.md).
