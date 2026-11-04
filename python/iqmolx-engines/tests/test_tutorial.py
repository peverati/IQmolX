# Copyright (C) 2026 Roberto Peverati
# Developed with the assistance of Claude (Anthropic).
#
# This file is part of IQMolX, a fork of IQmol (http://iqmol.org). It is free
# software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version. It is
# distributed WITHOUT ANY WARRANTY; see <http://www.gnu.org/licenses/>.

"""The in-app tutorial (resources/help/plugin-tutorial.html) shows code of the tutorial plugin.
Every block of code that is copied from a file is marked in the page with
<!-- file: PATH FIRST-LAST --> and must be exactly those lines of the file, so that the tutorial
cannot drift from the code that the tests run.  The images and links of the page must exist."""

import html
import os
import re

HERE = os.path.dirname(__file__)
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
PAGE = os.path.join(ROOT, "resources", "help", "plugin-tutorial.html")


def page():
    with open(PAGE, encoding="utf-8") as f:
        return f.read()


def test_every_marked_block_is_the_code_of_its_file():
    blocks = re.findall(r"<!-- file: (\S+) (\d+)-(\d+) -->\n<pre>(.*?)</pre>", page(), re.S)
    assert len(blocks) >= 7, "the tutorial shows the code of the plugin in blocks"
    for path, first, last, shown in blocks:
        with open(os.path.join(ROOT, path), encoding="utf-8") as f:
            wanted = "\n".join(f.read().split("\n")[int(first) - 1:int(last)])
        assert html.unescape(shown) == wanted, "%s lines %s-%s changed: run scripts/make-plugin-tutorial.py" % (path, first, last)


def test_the_whole_plugin_is_shown():
    """Every line of the engine class and its methods is in some block."""
    path = "python/examples/iqmolx-engine-emt/src/iqmolx_engine_emt/__init__.py"
    covered = set()
    for p, first, last in re.findall(r"<!-- file: (\S+) (\d+)-(\d+) -->", page()):
        if p == path:
            covered.update(range(int(first), int(last) + 1))
    with open(os.path.join(ROOT, path), encoding="utf-8") as f:
        lines = f.read().split("\n")
    start = next(i for i, l in enumerate(lines, 1) if l.startswith("import importlib"))
    code = [i for i in range(start, len(lines) + 1) if lines[i - 1].strip()]
    missing = [i for i in code if i not in covered]
    assert not missing, "lines of the plugin that the tutorial does not show: %s" % missing


def test_images_and_links_exist_and_are_in_the_resources():
    text = page()
    helpdir = os.path.dirname(PAGE)
    qrc = open(os.path.join(ROOT, "resources", "IQmol.qrc"), encoding="utf-8").read()
    for image in set(re.findall(r'<img src="([^"]+)"', text)):
        assert os.path.exists(os.path.join(helpdir, image)), image
        assert "help/" + image in qrc, image + " is not in IQmol.qrc"
    assert "help/plugin-tutorial.html" in qrc
    for target in set(re.findall(r'href="([a-z-]+\.html)', text)):
        assert os.path.exists(os.path.join(helpdir, target)), target
    names = set(re.findall(r'<a name="([^"]+)">', text))
    for anchor in re.findall(r'href="#([^"]+)"', text):
        assert anchor in names, anchor


def test_the_help_pages_link_to_the_tutorial():
    for name in ("index.html", "python-engines.html"):
        with open(os.path.join(os.path.dirname(PAGE), name), encoding="utf-8") as f:
            assert "plugin-tutorial.html" in f.read(), name
