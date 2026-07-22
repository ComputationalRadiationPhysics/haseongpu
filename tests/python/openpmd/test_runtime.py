from pathlib import Path

import pytest

from pyInclude import _runtime
from pyInclude.openpmd import transport


def _write_runtime_config(runtime_dir, *, provider_dir="/provider/one"):
    path = runtime_dir / "python" / "pyInclude" / "_runtime" / "_config.py"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f'HASE_RUNTIME_DIR = {str(runtime_dir)!r}\n'
        'HASE_RUNTIME_VERSION = "2.1.0"\n'
        'HASE_USE_SYSTEM_OPENPMD = True\n'
        f'HASE_OPENPMD_PYTHON_PACKAGE_DIR = {provider_dir!r}\n',
        encoding="utf-8",
    )
    return path


def test_nonstandardRuntimeOwnsExecutableLibrariesAndMetadata(tmp_path, monkeypatch):
    runtime_dir = tmp_path / "runtime-with-a-nonstandard-name"
    _write_runtime_config(runtime_dir)
    monkeypatch.setenv("HASE_RUNTIME_DIR", str(runtime_dir))

    assert _runtime.runtime_root() == runtime_dir
    assert _runtime.runtime_config().HASE_OPENPMD_PYTHON_PACKAGE_DIR == "/provider/one"
    assert next(_runtime.runtime_executable_candidates(("calcPhiASE",))).parent == runtime_dir
    assert next(_runtime.runtime_library_candidates(("libprobe.so",))).parent == runtime_dir / "lib"


def test_runtimeMetadataIsReloadedAfterNativeReconfigure(tmp_path, monkeypatch):
    runtime_dir = tmp_path / "build"
    _write_runtime_config(runtime_dir, provider_dir="/provider/one")
    monkeypatch.setenv("HASE_RUNTIME_DIR", str(runtime_dir))

    assert _runtime.runtime_config().HASE_OPENPMD_PYTHON_PACKAGE_DIR == "/provider/one"
    _write_runtime_config(runtime_dir, provider_dir="/provider/two")
    assert _runtime.runtime_config().HASE_OPENPMD_PYTHON_PACKAGE_DIR == "/provider/two"


def test_missingRuntimeMetadataIsNotReportedAsProviderIncompatibility(tmp_path, monkeypatch):
    runtime_dir = tmp_path / "unconfigured-runtime"
    runtime_dir.mkdir()
    monkeypatch.setenv("HASE_RUNTIME_DIR", str(runtime_dir))

    with pytest.raises(RuntimeError, match="runtime metadata is missing") as error:
        transport._runtime_config()

    assert "incompatible" not in str(error.value).lower()
    assert str(_runtime.runtime_metadata_path(runtime_dir)) in str(error.value)


def test_sourceFrontendUsesProviderMetadataFromSelectedRuntime(tmp_path, monkeypatch):
    runtime_dir = tmp_path / "build-tuned-under-an-arbitrary-name"
    provider_dir = tmp_path / "provider" / "site-packages"
    (provider_dir / "openpmd_api").mkdir(parents=True)
    _write_runtime_config(runtime_dir, provider_dir=str(provider_dir))
    monkeypatch.setenv("HASE_RUNTIME_DIR", str(runtime_dir))
    monkeypatch.delenv("HASE_OPENPMD_PYTHONPATH", raising=False)
    monkeypatch.delenv("HASE_OPENPMD_PYTHON_PACKAGE_DIR", raising=False)

    candidates = list(transport._candidate_python_paths(runtime_dir / "calcPhiASE"))

    assert provider_dir in candidates
