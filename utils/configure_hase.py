#!/usr/bin/env python3
"""Source-tree wrapper for the installed hase-configure command."""

from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path


def _load_source_tree_main():
    """Load pyInclude.config without requiring the compiled extension first."""
    repo_root = Path(__file__).resolve().parents[1]
    package = sys.modules.get("pyInclude")
    if package is None:
        package = types.ModuleType("pyInclude")
        package.__path__ = [str(repo_root / "pyInclude")]
        package.__package__ = "pyInclude"
        sys.modules["pyInclude"] = package

    spec = importlib.util.spec_from_file_location(
        "pyInclude.config",
        repo_root / "pyInclude" / "config.py",
        submodule_search_locations=None,
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["pyInclude.config"] = module
    spec.loader.exec_module(module)
    return module.main


if __name__ == "__main__":
    raise SystemExit(_load_source_tree_main()())
