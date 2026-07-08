import importlib.util
from pathlib import Path
from types import SimpleNamespace


def _load_preflight_module():
    root = Path(__file__).resolve().parents[3]
    path = root / "utils" / "check_openpmd_compatibility.py"
    spec = importlib.util.spec_from_file_location("check_openpmd_compatibility", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_pythonBackendCheckAcceptsDefaultAdiosSstProvider():
    preflight = _load_preflight_module()
    errors = []
    info = {
        "variants": {"adios2": True},
        "file_extensions": ["bp", "sst"],
    }

    preflight._check_python_backend(info, "adios-sst", errors)

    assert errors == []


def test_pythonBackendCheckRejectsMissingHdf5Support():
    preflight = _load_preflight_module()
    errors = []
    info = {
        "variants": {"adios2": True, "hdf5": False},
        "file_extensions": ["bp", "sst"],
    }

    preflight._check_python_backend(info, "hdf5", errors)

    assert "HDF5 support" in "\n".join(errors)


def test_backendSummaryMarksUnconfirmedSstSupport():
    preflight = _load_preflight_module()
    info = {
        "openPMD_VERSION": "0.17.0",
        "openPMD_HAVE_ADIOS2": "TRUE",
        "openPMD_HAVE_HDF5": "TRUE",
        "ADIOS2_HAVE_SST": "<undefined>",
    }

    assert preflight._cmake_supported_backends(info) == [
        "adios",
        "adios-sst (SST unconfirmed)",
        "hdf5",
    ]


def test_pythonBackendSummaryReportsSupportedRuntimeBackends():
    preflight = _load_preflight_module()
    info = {
        "variants": {"adios2": True, "hdf5": True},
        "file_extensions": ["bp", "h5", "sst"],
    }

    assert preflight._python_supported_backends(info) == [
        "adios",
        "adios-sst",
        "hdf5",
    ]


def test_defaultCmakeGeneratorPrefersNinja(monkeypatch, tmp_path):
    preflight = _load_preflight_module()
    ninja = tmp_path / "ninja"
    ninja.touch()
    ninja.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    monkeypatch.delenv("CMAKE_GENERATOR", raising=False)

    assert preflight._default_cmake_generator() == "Ninja"


def test_defaultCmakeGeneratorRespectsCmakeGeneratorEnv(monkeypatch, tmp_path):
    preflight = _load_preflight_module()
    ninja = tmp_path / "ninja"
    ninja.touch()
    ninja.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    monkeypatch.setenv("CMAKE_GENERATOR", "Unix Makefiles")

    assert preflight._default_cmake_generator() is None


def test_cmakeProbeReportsMissingCmakeExecutable():
    preflight = _load_preflight_module()
    errors = []
    info = preflight._cmake_probe(
        SimpleNamespace(
            cmake="/definitely/missing/cmake",
            cmake_generator=None,
            cmake_prefix_path=None,
            openpmd_dir=None,
        ),
        errors,
    )

    assert "was not found" in "\n".join(errors)
    assert "/definitely/missing/cmake" in info["command"]
