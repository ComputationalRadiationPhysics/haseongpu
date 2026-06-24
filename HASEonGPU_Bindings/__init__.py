# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

from importlib.machinery import PathFinder
from pathlib import Path
import sys

_package_dir = Path(__file__).resolve().parent
_repo_root = _package_dir.parent


def _build_binding_dirs():
    legacy_dir = _repo_root / "build" / "python" / "HASEonGPU_Bindings"
    if legacy_dir.is_dir():
        yield legacy_dir

    build_root = _repo_root / "build"
    if not build_root.is_dir():
        return

    for candidate in sorted(build_root.glob("*/python/HASEonGPU_Bindings")):
        if candidate.is_dir() and candidate != legacy_dir:
            yield candidate


def _non_repo_sys_path():
    filtered = []
    for entry in sys.path:
        entry_path = Path(entry or ".").resolve()
        if entry_path == _repo_root:
            continue
        filtered.append(entry)
    return filtered


def _installed_binding_dirs():
    spec = PathFinder.find_spec(__name__, _non_repo_sys_path())
    if spec is None:
        return

    for location in spec.submodule_search_locations or ():
        candidate = Path(location).resolve()
        if candidate != _package_dir:
            yield candidate


def _ordered_package_paths():
    seen = set()
    for candidate in (*_build_binding_dirs(), _package_dir, *_installed_binding_dirs()):
        normalized = str(candidate.resolve())
        if normalized in seen:
            continue
        seen.add(normalized)
        yield normalized


__path__[:] = list(_ordered_package_paths())

from .HASEonGPU import *
