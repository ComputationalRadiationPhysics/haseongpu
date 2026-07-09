import importlib.util
from pathlib import Path

import pytest


def _source_backends_module():
    root = Path(__file__).resolve().parents[3]
    path = root / "pyInclude" / "openpmd" / "backends.py"
    spec = importlib.util.spec_from_file_location("hase_source_openpmd_backends", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


backends = _source_backends_module()


class FakeProbeFunction:
    def __init__(self, callback):
        self.callback = callback
        self.argtypes = None
        self.restype = None

    def __call__(self, *args):
        return self.callback(*args)


class FakeProbeLibrary:
    def __init__(self, names):
        self.names = tuple(name.encode("utf-8") for name in names)
        self.haseOpenPmdBackendCount = FakeProbeFunction(lambda: len(self.names))
        self.haseOpenPmdBackendName = FakeProbeFunction(
            lambda index: self.names[index] if index < len(self.names) else None
        )


def test_openPmdBackendsLoadsCompiledProbe(monkeypatch, tmp_path):
    probe_path = tmp_path / "libHaseOpenPmdBackendProbe.so"
    monkeypatch.setattr(
        backends,
        "_load_probe_library",
        lambda extra_dirs=(): (FakeProbeLibrary(("hdf5", "adios", "unsupported")), probe_path),
    )

    assert backends._load_backend_names() == (("adios", "hdf5"), probe_path)


def test_openPmdBackendsCachesNames(monkeypatch):
    calls = []

    def fake_load_backend_names(extra_dirs=()):
        calls.append(tuple(extra_dirs))
        return ("adios",), "probe"

    monkeypatch.setattr(backends, "_load_backend_names", fake_load_backend_names)
    monkeypatch.setattr(backends.OpenPmdBackends, "_known", None)
    monkeypatch.setattr(backends.OpenPmdBackends, "_probe_path", None)

    assert backends.OpenPmdBackends.all() == ["adios"]
    assert backends.OpenPmdBackends.known() == ["adios"]
    assert calls == [()]


def test_openPmdBackendsFailsWhenProbeReportsNoKnownBackends(monkeypatch, tmp_path):
    probe_path = tmp_path / "libHaseOpenPmdBackendProbe.so"
    monkeypatch.setattr(
        backends,
        "_load_probe_library",
        lambda extra_dirs=(): (FakeProbeLibrary(("unsupported",)), probe_path),
    )

    with pytest.raises(RuntimeError, match="did not report any supported backends"):
        backends._load_backend_names()
