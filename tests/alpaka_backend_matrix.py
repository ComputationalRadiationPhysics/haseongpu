# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import os


RUNTIME_TEST_BACKEND_OVERRIDE = "HASE_TEST_ALPAKA_BACKEND"


def _compiled_backends():
    from HASEonGPU import AlpakaBackends

    return AlpakaBackends.all()


def alpaka_runtime_backend(backends=None, requested=None):
    available = list(_compiled_backends() if backends is None else backends)
    selection = (
        os.environ.get(RUNTIME_TEST_BACKEND_OVERRIDE, "").strip()
        if requested is None
        else requested.strip()
    )

    if selection == "cuda":
        matches = [
            backend
            for backend in available
            if "gpu" in backend.lower() and "cuda" in backend.lower()
        ]
        if len(matches) != 1:
            raise RuntimeError(
                "Expected exactly one CUDA GPU Alpaka backend, "
                f"found {matches}; all compiled backends: {available}"
            )
        return matches[0]

    if selection:
        if selection not in available:
            raise RuntimeError(
                f"Requested Alpaka backend '{selection}' is unavailable; "
                f"compiled backends: {available}"
            )
        return selection

    for preferred in ("Host_Cpu_CpuOmpBlocks", "Host_Cpu_CpuSerial"):
        if preferred in available:
            return preferred
    for backend in available:
        if "Cpu" in backend:
            return backend
    if available:
        return available[0]
    raise RuntimeError("no Alpaka backend is available in this build")


if __name__ == "__main__":
    print(alpaka_runtime_backend())
