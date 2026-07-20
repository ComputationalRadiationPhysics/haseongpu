# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Query Alpaka backend names compiled into the local HASEonGPU build."""

import ctypes
import ctypes.util
import sys
from pathlib import Path


def _libraryNames():
    if sys.platform == "win32":
        return ("HaseAlpakaBackendNames.dll",)
    if sys.platform == "darwin":
        return ("libHaseAlpakaBackendNames.dylib",)
    return ("libHaseAlpakaBackendNames.so",)


def _runtimePackageDirs():
    from importlib.util import find_spec

    spec = find_spec("pyInclude._runtime")
    if spec is None or spec.submodule_search_locations is None:
        return ()
    return tuple(Path(path) for path in spec.submodule_search_locations)


def _candidatePaths():
    moduleDir = Path(__file__).resolve().parent
    seen = set()

    def yieldPath(path):
        normalized = str(path.resolve())
        if normalized in seen:
            return
        seen.add(normalized)
        return path

    for packageDir in _runtimePackageDirs():
        for name in _libraryNames():
            candidate = yieldPath(packageDir / name)
            if candidate is not None:
                yield candidate

    for name in _libraryNames():
        candidate = yieldPath(moduleDir / name)
        if candidate is not None:
            yield candidate

    for parent in moduleDir.parents:
        for name in _libraryNames():
            candidate = yieldPath(parent / "build" / "python" / "pyInclude" / "_runtime" / name)
            if candidate is not None:
                yield candidate
        buildRoot = parent / "build"
        if not buildRoot.is_dir():
            continue
        for bindingDir in sorted(buildRoot.glob("*/python/pyInclude/_runtime")):
            for name in _libraryNames():
                candidate = yieldPath(bindingDir / name)
                if candidate is not None:
                    yield candidate
        for buildDir in sorted(buildRoot.glob("cp*")):
            for name in _libraryNames():
                candidate = yieldPath(buildDir / name)
                if candidate is not None:
                    yield candidate


def _loadLibrary():
    checked = []
    for path in _candidatePaths():
        checked.append(str(path))
        if path.is_file():
            return ctypes.CDLL(str(path))

    found = ctypes.util.find_library("HaseAlpakaBackendNames")
    if found is not None:
        return ctypes.CDLL(found)

    searched = "\n".join(f"  - {path}" for path in checked)
    raise RuntimeError(
        "Could not find the CMake-built Alpaka backend-name library. "
        "Build the HaseAlpakaBackendNames target first. Searched:\n" + searched
    )


def _loadBackendNames():
    lib = _loadLibrary()
    lib.haseAlpakaBackendCount.argtypes = ()
    lib.haseAlpakaBackendCount.restype = ctypes.c_size_t
    lib.haseAlpakaBackendName.argtypes = (ctypes.c_size_t,)
    lib.haseAlpakaBackendName.restype = ctypes.c_char_p

    names = []
    for index in range(lib.haseAlpakaBackendCount()):
        value = lib.haseAlpakaBackendName(index)
        if value is None:
            continue
        names.append(value.decode("utf-8"))
    return tuple(names)


class AlpakaBackends:
    """Namespace of backend strings accepted by ``PhiASE`` and ``calcPhiASE``.

    Backend availability depends on how the C++ extension was compiled. Use
    ``AlpakaBackends.all()`` to inspect the strings for this environment.
    """

    _known = None

    @classmethod
    def all(cls):
        """Return all backend names reported by the compiled helper library."""
        if cls._known is None:
            cls._known = _loadBackendNames()
            for backendName in cls._known:
                if backendName.isidentifier():
                    setattr(cls, backendName, backendName)
        return list(cls._known)

    @classmethod
    def known(cls):
        """Alias for ``all()`` kept for discoverability."""
        return cls.all()
