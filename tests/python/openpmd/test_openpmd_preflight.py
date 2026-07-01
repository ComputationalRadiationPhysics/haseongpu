import importlib.util
from pathlib import Path


def _load_preflight_module():
    root = Path(__file__).resolve().parents[3]
    path = root / "utils" / "check_openpmd_compatibility.py"
    spec = importlib.util.spec_from_file_location("check_openpmd_compatibility", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_python_backend_check_accepts_default_adios_sst_provider():
    preflight = _load_preflight_module()
    errors = []
    info = {
        "variants": {"adios2": True},
        "file_extensions": ["bp", "sst"],
    }

    preflight._check_python_backend(info, "adios-sst", errors)

    assert errors == []


def test_python_backend_check_rejects_missing_hdf5_support():
    preflight = _load_preflight_module()
    errors = []
    info = {
        "variants": {"adios2": True, "hdf5": False},
        "file_extensions": ["bp", "sst"],
    }

    preflight._check_python_backend(info, "hdf5", errors)

    assert "HDF5 support" in "\n".join(errors)


def test_backend_summary_marks_unconfirmed_sst_support():
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


def test_python_backend_summary_reports_supported_runtime_backends():
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
