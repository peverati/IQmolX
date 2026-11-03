import functools
import http.server
import os
import subprocess
import threading

import pytest

from iqmolx_engines import download
from iqmolx_engines.base import Cancelled, Context, EngineError


@pytest.fixture
def server(tmp_path):
    """A web server on this computer that serves tmp_path/www and counts the requests."""
    www = tmp_path / "www"
    www.mkdir()
    hits = []

    class Handler(http.server.SimpleHTTPRequestHandler):
        def do_GET(self):
            hits.append(self.path)
            super().do_GET()

        def log_message(self, *args):
            pass

    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0),
                                            functools.partial(Handler, directory=str(www)))
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield "http://127.0.0.1:%d" % httpd.server_address[1], www, hits
    httpd.shutdown()


def make_file(www, name="model.bin", size=3 * (1 << 20) + 123):
    data = os.urandom(size)
    (www / name).write_bytes(data)
    return str(www / name), size


class Recorder(Context):
    def __init__(self, cancel_after=None):
        self.events = []
        self.cancel_after = cancel_after
        super().__init__(progress=lambda f, m: self.events.append((f, m)),
                         cancelled=lambda: self.cancel_after is not None
                         and len(self.events) > self.cancel_after)


def test_blob_hash_is_the_one_git_computes(tmp_path):
    f = tmp_path / "x.bin"
    f.write_bytes(os.urandom(5000))
    git = subprocess.run(["git", "hash-object", str(f)], capture_output=True, text=True)
    if git.returncode != 0:
        pytest.skip("git is not available")
    assert download.git_blob_sha1(str(f)) == git.stdout.strip()


def test_download_with_progress_and_checks(server, tmp_path):
    url, www, hits = server
    src, size = make_file(www)
    dest = str(tmp_path / "cache" / "sub" / "model.bin")
    ctx = Recorder()
    got = download.fetch(url + "/model.bin", dest, ctx, label="the test model", size=size,
                         blob_sha1=download.git_blob_sha1(src))
    assert got == dest and open(dest, "rb").read() == open(src, "rb").read()
    fractions = [f for f, m in ctx.events if f is not None]
    assert fractions[0] == 0.0 and fractions[-1] == pytest.approx(1.0) and fractions == sorted(fractions)
    assert "the test model" in ctx.events[0][1]
    assert not os.path.exists(dest + ".part")


def test_a_matching_file_is_not_downloaded_again(server, tmp_path):
    url, www, hits = server
    src, size = make_file(www)
    dest = str(tmp_path / "model.bin")
    sha = download.git_blob_sha1(src)
    download.fetch(url + "/model.bin", dest, Context(), size=size, blob_sha1=sha)
    assert len(hits) == 1
    download.fetch(url + "/model.bin", dest, Context(), size=size, blob_sha1=sha)
    assert len(hits) == 1                                    # served from the cache


def test_a_wrong_checksum_or_size_is_refused_and_leaves_nothing(server, tmp_path):
    url, www, hits = server
    src, size = make_file(www)
    for kwargs in ({"blob_sha1": "0" * 40}, {"size": size + 1}):
        dest = str(tmp_path / "bad.bin")
        with pytest.raises(EngineError) as e:
            download.fetch(url + "/model.bin", dest, Context(), label="the model", **kwargs)
        assert e.value.code == "engine_error" and "not the expected one" in e.value.message
        assert not os.path.exists(dest) and not os.path.exists(dest + ".part")


def test_a_corrupt_cached_file_is_replaced(server, tmp_path):
    url, www, hits = server
    src, size = make_file(www)
    dest = str(tmp_path / "model.bin")
    open(dest, "wb").write(b"truncated")
    download.fetch(url + "/model.bin", dest, Context(), size=size, blob_sha1=download.git_blob_sha1(src))
    assert os.path.getsize(dest) == size and len(hits) == 1


def test_cancel_stops_the_download_and_removes_the_partial_file(server, tmp_path):
    url, www, hits = server
    make_file(www, size=8 * (1 << 20))
    dest = str(tmp_path / "model.bin")
    with pytest.raises(Cancelled):
        download.fetch(url + "/model.bin", dest, Recorder(cancel_after=2))
    assert not os.path.exists(dest) and not os.path.exists(dest + ".part")


def test_an_unreachable_server_is_reported(tmp_path):
    dest = str(tmp_path / "model.bin")
    with pytest.raises(EngineError) as e:
        download.fetch("http://127.0.0.1:9/model.bin", dest, Context(), label="the model")
    assert "Cannot download the model" in e.value.message and "127.0.0.1:9" in e.value.detail
    assert not os.path.exists(dest + ".part")


def test_cache_folder(monkeypatch, tmp_path):
    monkeypatch.setenv("IQMOLX_ENGINES_CACHE", str(tmp_path / "mine"))
    assert download.cache_dir() == str(tmp_path / "mine")
    monkeypatch.delenv("IQMOLX_ENGINES_CACHE")
    assert download.cache_dir().endswith(("Caches/IQMolX/engines", "iqmolx-engines"))
