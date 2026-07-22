# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Query Alpaka backend names compiled into the local HASEonGPU build."""

import ctypes
import ctypes.util
import sys

def _libraryNames():
    if sys.platform == "win32":
        return ("HaseAlpakaBackendNames.dll",)
    if sys.platform == "darwin":
        return ("libHaseAlpakaBackendNames.dylib",)
    return ("libHaseAlpakaBackendNames.so",)


def _candidatePaths():
    from pyInclude._runtime import runtime_library_candidates

    yield from runtime_library_candidates(_libraryNames())


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
