# Copyright (C) 2026 Roberto Peverati
# Developed with the assistance of Claude (Anthropic).
#
# This file is part of IQMolX, a fork of IQmol (http://iqmol.org). It is free
# software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version. It is
# distributed WITHOUT ANY WARRANTY; see <http://www.gnu.org/licenses/>.
"""Download of model files for engines, with progress, cancel, an atomic write and a check of
what arrived.  Engines whose weights are plain files (not fetched by their own library) use this;
the files are kept in a cache folder of the user, not in the program or the Python environment."""

import hashlib
import os
import sys
import urllib.error
import urllib.request

from .base import EngineError

CHUNK = 1 << 20


def cache_dir():
    """Where downloaded model files are kept: $IQMOLX_ENGINES_CACHE, else the cache folder of the
    system (~/Library/Caches/IQMolX/engines on a Mac, %LOCALAPPDATA%\\IQMolX\\engines on Windows)."""
    given = os.environ.get("IQMOLX_ENGINES_CACHE")
    if given:
        return given
    if sys.platform == "darwin":
        return os.path.expanduser("~/Library/Caches/IQMolX/engines")
    if sys.platform == "win32":
        base = os.environ.get("LOCALAPPDATA") or os.path.expanduser("~\\AppData\\Local")
        return os.path.join(base, "IQMolX", "engines")
    base = os.environ.get("XDG_CACHE_HOME") or os.path.expanduser("~/.cache")
    return os.path.join(base, "iqmolx-engines")


def git_blob_sha1(path):
    """The hash git uses for a file ("blob <size>\\0<content>"); GitHub's API lists it for every file,
    so a download from a repository can be checked without downloading it first."""
    sha = hashlib.sha1()
    sha.update(b"blob %d\0" % os.path.getsize(path))
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(CHUNK), b""):
            sha.update(chunk)
    return sha.hexdigest()


def _matches(path, size, blob_sha1):
    if not os.path.isfile(path):
        return False
    if size is not None and os.path.getsize(path) != size:
        return False
    return blob_sha1 is None or git_blob_sha1(path) == blob_sha1


def fetch(url, dest, ctx, label="the model", size=None, blob_sha1=None):
    """Make sure `dest` holds the file at `url`, downloading it if needed, and return its path.
    `size` and `blob_sha1` (the git blob hash) are checked when given; a file that does not match is
    deleted and the error says so.  A file that is already there and matches is not downloaded again."""
    if _matches(dest, size, blob_sha1):
        return dest
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    part = dest + ".part"
    request = urllib.request.Request(url, headers={"User-Agent": "iqmolx-engines"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response, open(part, "wb") as out:
            total = int(response.headers.get("Content-Length") or size or 0)
            done = 0
            ctx.progress(0.0 if total else None, "Downloading %s" % label)
            while True:
                ctx.check_cancel()
                chunk = response.read(CHUNK)
                if not chunk:
                    break
                out.write(chunk)
                done += len(chunk)
                ctx.progress(done / total if total else None,
                             "Downloading %s (%.1f of %.1f MB)" % (label, done / 1e6, total / 1e6)
                             if total else "Downloading %s (%.1f MB)" % (label, done / 1e6))
    except urllib.error.URLError as exc:
        _remove(part)
        raise EngineError("engine_error", "Cannot download %s: %s" % (label, exc.reason),
                          detail="Source: %s\nCheck the internet connection." % url) from None
    except OSError as exc:
        _remove(part)
        raise EngineError("engine_error", "Cannot download %s: %s" % (label, exc),
                          detail="Source: %s" % url) from None
    except BaseException:
        _remove(part)          # cancelled, or anything else: no half file is left behind
        raise
    if not _matches(part, size, blob_sha1):
        _remove(part)
        raise EngineError("engine_error", "The downloaded file of %s is not the expected one" % label,
                          detail="Source: %s (size or checksum differs from what is pinned)" % url)
    os.replace(part, dest)
    return dest


def _remove(path):
    try:
        os.remove(path)
    except OSError:
        pass
