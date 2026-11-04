# Copyright (C) 2026 Roberto Peverati
# Developed with the assistance of Claude (Anthropic).
#
# This file is part of IQMolX, a fork of IQmol (http://iqmol.org). It is free
# software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version. It is
# distributed WITHOUT ANY WARRANTY; see <http://www.gnu.org/licenses/>.

"""Writes resources/help/plugin-tutorial.html, the in-app tutorial on writing an engine plugin.

The code of the example plugin (python/examples/iqmolx-engine-emt) is copied from its files, with a
marker <!-- file: PATH FIRST-LAST --> before each block, and python/iqmolx-engines/tests/test_tutorial.py
checks that the blocks are still the code of the files.  When that test fails after a change of the
example, run this script again and read the page:

    python3 scripts/make-plugin-tutorial.py
"""
import html, os, re, sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
EXAMPLE = "python/examples/iqmolx-engine-emt"
INIT = EXAMPLE + "/src/iqmolx_engine_emt/__init__.py"
TEST = "python/iqmolx-engines/tests/test_plugin_emt.py"


def lines(path):
    with open(os.path.join(ROOT, path), encoding="utf-8") as f:
        return f.read().split("\n")


def block(path, first, last):
    """A marked, escaped <pre> with lines first..last (1-based, inclusive) of a file."""
    text = "\n".join(lines(path)[first - 1:last])
    return "<!-- file: %s %d-%d -->\n<pre>%s</pre>\n" % (path, first, last, html.escape(text, quote=False))


def find(path, pattern, after=0):
    for i, line in enumerate(lines(path), 1):
        if i > after and re.search(pattern, line):
            return i
    raise SystemExit("not found: " + pattern)


def pre(text):
    return "<pre>%s</pre>\n" % html.escape(text.strip("\n"), quote=False)


init = lambda a, b: block(INIT, a, b)
L = {
    "cls":      find(INIT, r"^class EmtEngine"),
    "requires": find(INIT, r"requirements = "),
    "options":  find(INIT, r"^    options = \("),
    "optend":   find(INIT, r"^    \)$"),
    "unavail":  find(INIT, r"@classmethod"),
    "loaddef":  find(INIT, r"def load\("),
    "compdef":  find(INIT, r"def compute\("),
    "end":      len(lines(INIT)) - 1 if lines(INIT)[-1] == "" else len(lines(INIT)),
}
L["loadend"] = L["compdef"] - 2
t0 = find(TEST, r"^def test_energy_and_forces_of_a_copper_dimer")
t1 = find(TEST, r"^def test_forces_are_minus", t0) - 3
py_end = find(EXAMPLE + "/pyproject.toml", r"^where = ")

page = """<html>
<head>
   <link href="style.css" rel="stylesheet" type="text/css" />
</head>

<body>

   <center> <h1> Writing an Engine Plugin: a Tutorial <br> </h1> </center><p></p>

   <p> A plugin adds a calculation method to the
   <a href="python-engines.html">Python Engines</a> window without changing or recompiling
   IQMolX: it is a small Python package, and the window builds itself from what the package
   says about its method.  This tutorial makes a complete plugin in about twenty minutes. It
   takes a calculator of the
   <a href="https://wiki.fysik.dtu.dk/ase/">Atomic Simulation Environment</a> (ASE), the
   effective medium theory (EMT), and turns it into an engine of its own with a name, a
   description, the right elements and settings.  The same steps work for any other
   method that you can call from Python.

   <p> You need the Python environment of IQMolX (see
   <a href="python-engines.html#setup">Setting Up Python</a>), a text editor and, for step 8,
   a terminal.  No programming beyond reading Python is needed to follow along.  The complete
   plugin is in the IQMolX sources, in <tt>python/examples/iqmolx-engine-emt</tt>, and the
   tests of IQMolX install and run it, so the code on this page is the code that works.  The
   reference for every detail is the developer guide <tt>PLUGINS.md</tt> in the
   <tt>python/iqmolx-engines</tt> folder.


   <br>
   <h3><a name="what"></a>1. What a Plugin Is</h3>

   <p> The calculations of the Python Engines window are done by a helper program that
   IQMolX starts in the Python environment.  An <i>engine</i> is a class in that helper.
   The class has two kinds of parts:
   <ul>
      <li> <b>attributes</b> that describe the method (its name, the elements it knows,
           whether it takes a charge, what it can calculate, its settings).  IQMolX reads
           them and builds its window from them, which is why it needs no code of its own
           for any engine; <br>
      <li> two <b>methods</b>, <tt>load</tt> (get the method ready) and <tt>compute</tt>
           (calculate the energy and forces of one molecule).
   </ul>
   A plugin is an ordinary Python package that contains such a class and announces it under
   the name <tt>iqmolx_engines</tt>.  The helper looks for these announcements every time the
   window asks for the engines, so a package that is installed while IQMolX is open is found
   with the <b>Find Engines</b> button.


   <br>
   <h3><a name="try"></a>2. First Try the Generic ASE Engine</h3>

   <p> Before writing code, see what you get for free.  The window has an engine that
   runs any ASE calculator: choose <i>ASE calculator</i> in the Engine list, press its
   <b>Install ase</b> button if it asks for it, and type
   <tt>ase.calculators.emt:EMT</tt> in <i>Calculator class</i>.

   <p><center><img src="images/PluginEngineAse.png" width="470" /></center></p>

   <p> Open <tt>cu4.xyz</tt> (from the example folder) with File&rarr;Open and run a geometry
   optimization: it works.  But the window knows nothing about EMT.  It says that all elements
   are possible (try iron: EMT has no parameters for it, and ASE answers with an error), the
   settings are a line of text and a JSON field, and there is no description, license or
   citation.  A plugin of its own fixes all of this.


   <br>
   <h3><a name="folder"></a>3. The Package</h3>

   <p> Make a folder <tt>iqmolx-engine-emt</tt> with this content (the names of the files
   matter, the name of the folder does not):

""" + pre("""
iqmolx-engine-emt/
    pyproject.toml
    src/
        iqmolx_engine_emt/
            __init__.py
""") + """
   <p> <tt>pyproject.toml</tt> describes the package and, most important, says where the
   engine is:

""" + block(EXAMPLE + "/pyproject.toml", 1, py_end) + """
   <ul>
      <li> <tt>[project.entry-points."iqmolx_engines"]</tt> is what makes the helper find the
           plugin.  The part after the <tt>=</tt> is <tt>module:Class</tt>.  The name on the
           left is free. <br>
      <li> <tt>dependencies</tt> lists what the plugin adds to the environment, here ASE.
           <b>Do not list <tt>iqmolx-engines</tt></b>: it is already in the environment and
           is not on PyPI, so the installer would not find it. <br>
      <li> A plugin that needs <tt>torch</tt> or a model package lists it here in the same
           way.
   </ul>


   <br>
   <h3><a name="class"></a>4. The Engine Class</h3>

   <p> <tt>__init__.py</tt> starts with the imports and the class.  Everything in it is
   what the window shows or checks:

""" + init(find(INIT, r"^import importlib"), L["cls"] - 2) + """
""" + init(L["cls"], L["requires"]) + """
   <table border="0" cellpadding="4">
      <tr><td valign="top"><tt>name</tt></td>
          <td> The identifier of the engine in the protocol and in the stored settings.  It
               must be different from the built-in engines (<tt>demo</tt>, <tt>ase</tt>,
               <tt>torchani</tt>, <tt>pyscf</tt>, ...): those win. </td></tr>
      <tr><td valign="top"><tt>title</tt>, <tt>description</tt></td>
          <td> The name in the Engine list and the text below it. </td></tr>
      <tr><td valign="top"><tt>license</tt>, <tt>citation</tt></td>
          <td> Shown in the window.  State the license of the method and of its parameters
               honestly; the window is where students look. </td></tr>
      <tr><td valign="top"><tt>elements</tt></td>
          <td> The atomic numbers that the method can do, or <tt>None</tt> for all.  The helper
               refuses other elements before your code runs, and the window shows the list. </td></tr>
      <tr><td valign="top"><tt>charge</tt>, <tt>multiplicity</tt></td>
          <td> Whether the method accepts a charged or open-shell molecule.  With
               <tt>False</tt> the helper refuses them. </td></tr>
      <tr><td valign="top"><tt>properties</tt></td>
          <td> What <tt>compute</tt> can return.  <tt>energy</tt> and <tt>forces</tt> are
               required.  Without a <tt>hessian</tt> IQMolX makes the vibrational frequencies from
               differences of forces by itself, as it does here. </td></tr>
      <tr><td valign="top"><tt>requirements</tt></td>
          <td> The pip requirements of the engine.  If they are missing, the engine is still
               listed and the window shows an <b>Install</b> button for exactly these. </td></tr>
   </table>


   <br>
   <h3><a name="options"></a>5. Settings</h3>

   <p> <tt>options</tt> describes the settings that appear in the window.  Each is made with
   <tt>option(name, title, type, default, ...)</tt>; the types are <tt>bool</tt>,
   <tt>int</tt>, <tt>float</tt>, <tt>string</tt> and <tt>choice</tt>.

""" + init(L["options"], L["optend"]) + """
   <p> The <i>scope</i> says when a setting is used.  A <b>load</b> setting (the default) is
   given to <tt>load</tt>, so changing it loads the method again; this is right for what is
   expensive to change, such as the model.  A <b>compute</b> setting is sent with every
   calculation to <tt>compute</tt>, so changing it costs nothing.  Settings can also have a
   <tt>min</tt> and <tt>max</tt>, a <tt>choices</tt> list, <tt>required</tt>,
   <tt>placeholder</tt> and <tt>show_if</tt> (shown only for some values of another setting);
   see the developer guide.

   <p><center><img src="images/PluginEngineEmt.png" width="470" /></center></p>

   <p> This is how the two settings appear in the window once the plugin is installed
   (step 7).


   <br>
   <h3><a name="available"></a>6. What Happens Without ASE</h3>

""" + init(L["unavail"], L["loaddef"] - 2) + """
   <p> A method needs its Python packages, and the window must be able to list the engine even
   when they are missing, so that it can offer to install them.  That is why ASE is not imported
   at the top of the file, and why <tt>unavailable_reason</tt> says what is missing without
   importing it.  An empty text means that the engine can run.


   <br>
   <h3><a name="work"></a>7. The Work: load and compute</h3>

""" + init(L["loaddef"], L["loadend"]) + """
   <p> <tt>load</tt> gets the method ready and is called once, when the engine is chosen or a
   <i>load</i> setting changes.  A neural network would load its weights here, and could
   report how far it is with <tt>ctx.progress(fraction, message)</tt>.  What it returns is
   shown in the status line.

""" + init(L["compdef"], L["end"]) + """
   <p> <tt>compute</tt> is called for every energy and force that IQMolX needs, and in a
   geometry optimization that is many times.  It receives the atomic numbers, the
   coordinates in <b>&Aring;ngstrom</b> (a numpy array), the charge and spin multiplicity, the
   list of what is wanted and the <i>compute</i> settings, and returns a dictionary with the
   energy in <b>eV</b> and the forces in <b>eV/&Aring;ngstrom</b> (minus the gradient).  These are
   the units of ASE, so nothing is converted here; for a method in other units the conversion
   belongs in this method.

   <p> A calculation that cannot be done is reported with
   <tt>raise EngineError(code, message, detail=...)</tt>.  The message appears in the window
   and the detail is kept for the log.  Any other Python error also becomes an error message
   instead of stopping IQMolX.  A long calculation should call <tt>ctx.check_cancel()</tt>
   now and then, so that the Cancel button works.


   <br>
   <h3><a name="install"></a>8. Install It and Use It</h3>

   <p> In the Python Engines window press <b>Install Plugin...</b>.

   <p><center><img src="images/PluginButtons.png" width="470" /></center></p>

   <p> A small window asks what to install: the name of a package (a plugin that someone
   published), a Git address, or, with <b>Folder...</b> and <b>File...</b>, a folder or a
   wheel on your computer.  Choose the folder <tt>iqmolx-engine-emt</tt> and press
   <b>Install</b>.

   <p><center><img src="images/PluginInstall.png" width="400" /></center></p>

   <p> The window runs the package installer of the environment, shows its output, and lists
   the engines again.  <i>EMT metals (ASE, example plugin)</i> is now in the Engine list
   with its description and settings.  Open <tt>cu4.xyz</tt>, choose the engine and run a
   geometry optimization.  Then move two atoms almost onto each other and run again: the
   window refuses with the message of your <tt>EngineError</tt>, and lowering <i>Closest
   atoms</i> to 0 lets the (meaningless) calculation through.

   <p> <b>Only install plugins from sources that you trust.</b>  Installing a package runs
   its code, and an engine runs on your computer with your files.

   <p> The same can be done in a terminal, which is also what you want while developing,
   because <tt>-e</tt> installs the folder in place, so that edits are picked up:

""" + pre("""
PY="$HOME/Library/Application Support/IQmol/python/bin/python"
uv pip install --python "$PY" -e ./iqmolx-engine-emt
""") + """
   <p> (<tt>uv</tt> is the installer that IQMolX uses; <tt>PY</tt> is the Python of the
   IQMolX environment.  If you chose your own Python in the window, use that one.)  Then press
   <b>Find Engines</b>: a running helper notices a plugin that was installed after it started.
   <b>The code of a plugin that was already loaded is not read again</b>: after you edit it,
   quit and start IQMolX.


   <br>
   <h3><a name="test"></a>9. Test It Without IQMolX</h3>

   <p> The helper can be asked what it finds, with nothing else running.  This lists your
   engine and its description as the window will see it:

""" + pre("""
"$PY" -m iqmolx_engines --describe | grep -A3 '"name": "emt"'
""") + """
   <p> The helper speaks one JSON message per line on its standard input and output (the
   <tt>PROTOCOL.md</tt> in the same folder describes every message), so a calculation can
   also be tried by hand.  Load the engine, calculate a copper dimer at 2.4 &Aring;, and then
   two atoms that almost touch:

""" + pre("""
printf '%s\\n' \\
  '{"id":1,"type":"load","engine":"emt"}' \\
  '{"id":2,"type":"compute","numbers":[29,29],"coordinates":[[0,0,0],[0,0,2.4]]}' \\
  '{"id":3,"type":"compute","numbers":[29,29],"coordinates":[[0,0,0],[0,0,0.05]]}' \\
  | "$PY" -m iqmolx_engines 2>/dev/null
""") + """
   <p> The answers are (the first line announces the helper):

""" + pre("""
{"event":"ready","protocol":1,"helper":"0.1.1","python":"3.12.7"}
{"id":1,"ok":true,"result":{"engine":"emt","info":{"calculator":"ase.calculators.emt.EMT"}}}
{"id":2,"ok":true,"result":{"energy":3.504405458151052,"forces":[[0.0,0.0,2.3571882814086518],[0.0,0.0,-2.3571882814086518]]}}
{"id":3,"ok":false,"error":{"code":"engine_error","message":"Two atoms are only 0.05 Angstrom apart; EMT is not meaningful for that","detail":"Move the atoms apart, or lower the closest distance allowed in the settings."}}
""") + """
   <p> For a real plugin, write tests with <tt>pytest</tt>.  The test file of this example,
   <tt>tests/test_plugin_emt.py</tt> of the helper, calls the engine directly and checks the
   energy against ASE, the forces against differences of the energy, the error, and the
   whole way through the installed plugin.  This is the first one:

""" + block(TEST, t0, t1) + """
   <p> Forces that are not minus the gradient of the energy are the commonest mistake of a
   new engine (a geometry optimization then wanders or never converges), so a test of
   the forces by finite differences, like the second test of that file, is worth copying.


   <br>
   <h3><a name="share"></a>10. Share It</h3>

   <ul>
      <li> Put the folder on GitHub (or send it as a zip); others install it with
           <b>Install Plugin...</b> from the folder, the wheel or the Git address
           (<tt>git+https://github.com/you/iqmolx-engine-emt</tt>). <br>
      <li> To publish it on PyPI, build it with <tt>python -m build</tt> and upload it
           with <tt>twine</tt>; users then only type the package name. <br>
      <li> <b>Weights and parameters are never bundled with IQMolX.</b>  A method that needs
           a model file downloads it the first time it is used, with a progress message, and
           says in <tt>license</tt> what applies to it.  Check the license before you
           redistribute anything.  The developer guide has a section on this, and a checklist
           for a plugin that is ready to be shared.
   </ul>


   <br>
   <h3><a name="next"></a>11. Where to Go from Here</h3>

   <ul>
      <li> <b>Another ASE calculator.</b>  Change <tt>load</tt> and the settings: for
           example <tt>TBLite</tt> of the <tt>tblite</tt> package (GFN2-xTB with ASE) takes
           a charge, so <tt>charge = True</tt> and the charge arrives in <tt>compute</tt>
           (the generic ASE engine passes charge and spin on to calculators that have such
           arguments, which is worth reading in its source, <tt>ase_adapter.py</tt>). <br>
      <li> <b>A machine learning model.</b>  Load the weights in <tt>load</tt>, with
           <tt>ctx.progress</tt> for a download, set <tt>requirements</tt> to the packages, and
           convert units in <tt>compute</tt> if needed.  The ANI and Egret engines of IQMolX
           (<tt>engines/torchani_engine.py</tt>, <tt>engines/egret_engine.py</tt>) are real
           examples. <br>
      <li> <b>More properties.</b>  Add <tt>dipole</tt> (in e&middot;&Aring;), <tt>charges</tt> or
           an analytic <tt>hessian</tt> (eV/&Aring;<sup>2</sup>, 3N&times;3N) to
           <tt>properties</tt> and return them from <tt>compute</tt> when they are in
           <tt>want</tt>. <br>
      <li> <b>The details.</b>  <tt>PLUGINS.md</tt> explains every attribute, the error
           codes, progress and cancel, and has a checklist.
   </ul>


   <br>
   <h3><a name="problems"></a>If Something Goes Wrong</h3>

   <table border="0" cellpadding="4">
      <tr><td valign="top"><i>The engine is not in the list</i></td>
          <td> Press <b>Find Engines</b>.  Then check that the plugin was installed into the
               Python shown in the window (a plugin in another Python is not seen), and run
               the <tt>--describe</tt> command of step 9 with the same Python: the helper
               prints a line on its error output for a plugin that it cannot use, for
               example <tt>cannot load plugin emt: ...</tt> (the package fails to import),
               <tt>is not an Engine subclass</tt>, or <tt>has no engine name</tt>. </td></tr>
      <tr><td valign="top"><i>Not available</i> in the list</td>
          <td> A package of <tt>requirements</tt> is missing: press the <b>Install</b>
               button that appears next to the list. </td></tr>
      <tr><td valign="top"><i>hidden by the built-in engine</i></td>
          <td> The <tt>name</tt> of the plugin is the name of an engine of IQMolX.  Choose
               another. </td></tr>
      <tr><td valign="top">A change to the code has no effect</td>
          <td> The helper keeps what it has already imported.  Quit and start IQMolX. </td></tr>
      <tr><td valign="top">The installer reports an error</td>
          <td> Read its output in the window: a misspelled name, a missing
               <tt>pyproject.toml</tt>, or a dependency that cannot be installed.  A plugin
               must not list <tt>iqmolx-engines</tt> as a dependency. </td></tr>
      <tr><td valign="top">Optimizations do not converge</td>
          <td> Check that the forces are minus the gradient of the energy, in eV/&Aring;ngstrom
               (test 2 of the example). </td></tr>
   </table>

</body>

</html>
"""
out = os.path.join(ROOT, "resources/help/plugin-tutorial.html")
with open(out, "w", encoding="utf-8") as f:
    f.write(page)
print("written", out, len(page))
