# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import importlib.util
import os
from pathlib import Path

BACKEND_PRIORITY = ("adios-sst", "adios", "hdf5")
ARTIFACT_NAME = "openpmd_backends.txt"


def _clean(values):
    selected = {str(value).strip().lower() for value in values if str(value).strip()}
    return [backend for backend in BACKEND_PRIORITY if backend in selected]


def _env_backends():
    return _clean(os.environ.get("HASE_OPENPMD_TEST_BACKENDS", "").split(","))


def _binding_package_dirs():
    spec = importlib.util.find_spec("HASEonGPU_Bindings")
    if spec is None or spec.submodule_search_locations is None:
        return ()
    return tuple(Path(path) for path in spec.submodule_search_locations)


def _build_artifact_candidates():
    configured = os.environ.get("HASE_OPENPMD_BACKEND_PROBE") or os.environ.get("HASE_CALCPHIASE")
    if configured:
        path = Path(configured).resolve()
        yield path.with_name(ARTIFACT_NAME)
        yield path.parent / "python" / "HASEonGPU_Bindings" / ARTIFACT_NAME

    root = Path(__file__).resolve().parents[1]
    yield root / "build" / "python" / "HASEonGPU_Bindings" / ARTIFACT_NAME
    yield root / "build" / "ci" / "python" / "HASEonGPU_Bindings" / ARTIFACT_NAME
    for path in sorted(root.glob("build/*/python/HASEonGPU_Bindings/" + ARTIFACT_NAME)):
        yield path


def _artifact_backends():
    candidates = [package_dir / ARTIFACT_NAME for package_dir in _binding_package_dirs()]
    candidates.extend(_build_artifact_candidates())
    for path in candidates:
        if path.is_file():
            backends = _clean(path.read_text(encoding="utf-8").splitlines())
            if not backends:
                raise RuntimeError(f"openPMD backend probe artifact is empty: {path}")
            return backends
    searched = "\n".join(f"  - {path}" for path in candidates)
    raise RuntimeError(
        "No openPMD backend probe artifact found. Build HaseOpenPmdBackendProbe "
        "or set HASE_OPENPMD_TEST_BACKENDS for an explicit manual matrix. Searched:\n" + searched
    )


def openpmd_test_backends():
    return _env_backends() or _artifact_backends()


def openpmd_runtime_backend():
    return openpmd_test_backends()[0]
