#!/usr/bin/env python3
"""Check that Python and CMake openPMD providers match HASEonGPU needs."""

from __future__ import annotations

import importlib.util
from pathlib import Path


def _load_preflight_module():
    path = Path(__file__).resolve().parents[1] / "pyInclude" / "openpmd" / "preflight.py"
    spec = importlib.util.spec_from_file_location("hase_openpmd_preflight", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_preflight = _load_preflight_module()

BACKEND_EXTENSIONS = _preflight.BACKEND_EXTENSIONS
HASE_SOURCE_BUILD_BACKENDS = _preflight.HASE_SOURCE_BUILD_BACKENDS
_check_cmake_backend = _preflight._check_cmake_backend
_check_python_backend = _preflight._check_python_backend
_check_versions = _preflight._check_versions
_cmake_bool = _preflight._cmake_bool
_cmake_probe = _preflight._cmake_probe
_cmake_supported_backends = _preflight._cmake_supported_backends
_default_cmake_generator = _preflight._default_cmake_generator
_format_backend_list = _preflight._format_backend_list
_parse_args = _preflight._parse_args
_print_backend_summary = _preflight._print_backend_summary
_python_info = _preflight._python_info
_python_supported_backends = _preflight._python_supported_backends
_version_tuple = _preflight._version_tuple
check_provider = _preflight.check_provider
main = _preflight.main
print_report = _preflight.print_report

__all__ = [
    "BACKEND_EXTENSIONS",
    "HASE_SOURCE_BUILD_BACKENDS",
    "_check_cmake_backend",
    "_check_python_backend",
    "_check_versions",
    "_cmake_bool",
    "_cmake_probe",
    "_cmake_supported_backends",
    "_default_cmake_generator",
    "_format_backend_list",
    "_parse_args",
    "_print_backend_summary",
    "_python_info",
    "_python_supported_backends",
    "_version_tuple",
    "check_provider",
    "main",
    "print_report",
]


if __name__ == "__main__":
    raise SystemExit(main())
