# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

from __future__ import annotations

import importlib
import sys
from pathlib import Path


def _resolve_sys_path_entry(entry):
    try:
        return (Path.cwd() if entry == "" else Path(entry)).resolve()
    except OSError:
        return Path(entry)


def ensure_hase_importable():
    """Prefer an installed HASEonGPU package, then fall back to this checkout."""
    source_root = Path(__file__).resolve().parents[1]
    original_path = list(sys.path)
    sys.path[:] = [
        entry
        for entry in original_path
        if _resolve_sys_path_entry(entry) != source_root
    ]
    try:
        importlib.import_module("HASEonGPU")
        return
    except ModuleNotFoundError as err:
        if err.name not in {
            "HASEonGPU",
            "HASEonGPU_Bindings",
            "HASEonGPU_Bindings.HASEonGPU",
        }:
            raise
    finally:
        sys.path[:] = original_path

    source_root_text = str(source_root)
    if source_root_text not in sys.path:
        sys.path.insert(0, source_root_text)
