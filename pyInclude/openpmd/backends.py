# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Query openPMD backends available to the local HASEonGPU runtime."""

from __future__ import annotations

import ctypes
import ctypes.util
import os
import sys
from pathlib import Path

BACKEND_PRIORITY = ("adios-sst", "adios", "hdf5")


def _library_names():
    if sys.platform == "win32":
        return ("HaseOpenPmdBackendProbe.dll",)
    if sys.platform == "darwin":
        return ("libHaseOpenPmdBackendProbe.dylib",)
    return ("libHaseOpenPmdBackendProbe.so",)


def _candidate_paths(extra_dirs=()):
    from pyInclude._runtime import runtime_library_candidates

    seen = set()

    def yield_path(path):
        normalized = str(path.resolve())
        if normalized in seen:
            return None
        seen.add(normalized)
        return path

    configured = os.environ.get("HASE_OPENPMD_BACKEND_PROBE_LIBRARY")
    if configured:
        candidate = yield_path(Path(configured))
        if candidate is not None:
            yield candidate

    for directory in extra_dirs:
        for name in _library_names():
            candidate = yield_path(Path(directory) / name)
            if candidate is not None:
                yield candidate

    for runtime_candidate in runtime_library_candidates(_library_names()):
        candidate = yield_path(runtime_candidate)
        if candidate is not None:
            yield candidate


def _load_probe_library(extra_dirs=()):
    checked = []
    for path in _candidate_paths(extra_dirs):
        checked.append(str(path))
        if path.is_file():
            return ctypes.CDLL(str(path)), path

    found = ctypes.util.find_library("HaseOpenPmdBackendProbe")
    if found is not None:
        return ctypes.CDLL(found), Path(found)

    searched = "\n".join(f"  - {path}" for path in checked)
    raise RuntimeError(
        "Could not find the CMake-built openPMD backend-probe library. "
        "Build the HaseOpenPmdBackendProbe target first. Searched:\n" + searched
    )


def _clean_backend_names(values):
    selected = {str(value).strip().lower() for value in values if str(value).strip()}
    return tuple(backend for backend in BACKEND_PRIORITY if backend in selected)


def _load_backend_names(extra_dirs=()):
    lib, path = _load_probe_library(extra_dirs)
    lib.haseOpenPmdBackendCount.argtypes = ()
    lib.haseOpenPmdBackendCount.restype = ctypes.c_size_t
    lib.haseOpenPmdBackendName.argtypes = (ctypes.c_size_t,)
    lib.haseOpenPmdBackendName.restype = ctypes.c_char_p

    names = []
    for index in range(lib.haseOpenPmdBackendCount()):
        value = lib.haseOpenPmdBackendName(index)
        if value is None:
            continue
        names.append(value.decode("utf-8"))

    backends = _clean_backend_names(names)
    if not backends:
        raise RuntimeError(f"openPMD backend-probe library did not report any supported backends: {path}")
    return backends, path


class OpenPmdBackends:
    """Namespace of openPMD backend strings accepted by HASEonGPU."""

    _known = None
    _probe_path = None

    @classmethod
    def all(cls):
        """Return all openPMD backends reported by the compiled probe library."""
        if cls._known is None:
            cls._known, cls._probe_path = _load_backend_names()
        return list(cls._known)

    @classmethod
    def known(cls):
        """Alias for ``all()`` kept for discoverability."""
        return cls.all()


for backend_name in BACKEND_PRIORITY:
    attr = backend_name.replace("-", "_")
    if attr.isidentifier():
        setattr(OpenPmdBackends, attr, backend_name)
