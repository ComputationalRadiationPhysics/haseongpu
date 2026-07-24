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


def _clear_hase_modules():
    for name in list(sys.modules):
        if name == "HASEonGPU" or name.startswith("HASEonGPU."):
            del sys.modules[name]
        elif name == "pyInclude" or name.startswith("pyInclude."):
            del sys.modules[name]


def ensure_hase_importable():
    """Prefer a compatible installed HASE package, then fall back to this checkout."""
    source_root = Path(__file__).resolve().parents[1]
    original_path = list(sys.path)
    sys.path[:] = [
        entry
        for entry in original_path
        if _resolve_sys_path_entry(entry) != source_root
    ]
    try:
        module = importlib.import_module("HASEonGPU")
        if hasattr(module, "VolumeTopology"):
            return
    except ModuleNotFoundError as err:
        if err.name != "HASEonGPU":
            raise
    finally:
        sys.path[:] = original_path

    _clear_hase_modules()
    source_root_text = str(source_root)
    if source_root_text not in sys.path:
        sys.path.insert(0, source_root_text)
