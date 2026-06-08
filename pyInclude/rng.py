# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""RNG seed helpers for Python-launched backend invocations."""

from __future__ import annotations

import threading

import numpy as np


_BACKEND_ROOT_SEED = int(np.random.SeedSequence().generate_state(1, dtype=np.uint32)[0])
_BACKEND_SEED_RNG = np.random.default_rng(_BACKEND_ROOT_SEED)
_BACKEND_SEED_LOCK = threading.Lock()


def backendRootRngSeed() -> int:
    """Return the non-deterministic seed used for this process-local seed stream."""
    return _BACKEND_ROOT_SEED


def defaultBackendRngSeed() -> int:
    """Return a fresh uint32 seed for one backend invocation."""
    with _BACKEND_SEED_LOCK:
        return int(_BACKEND_SEED_RNG.integers(0, np.iinfo(np.uint32).max + 1, dtype=np.uint32))
