# iqmolx-engine-emt

The example plugin of the tutorial *Writing an engine plugin* (in the IQMolX Help, and
`python/iqmolx-engines/PLUGINS.md` for the details): the effective medium theory calculator of
ASE as an engine of its own.  It is real code that the tests of `iqmolx-engines` install and run.

Install it from IQMolX with **Install Plugin...** in the Python Engines window (choose this
folder), or from a terminal into the environment of IQMolX:

    uv pip install --python "$HOME/Library/Application Support/IQmol/python/bin/python" ./iqmolx-engine-emt

Then press **Refresh** (the circular arrows next to the Engine list): *EMT metals (ASE, example plugin)* is in the Engine list.  `cu4.xyz`
is a starting geometry (File, Open, then run a geometry optimization).
