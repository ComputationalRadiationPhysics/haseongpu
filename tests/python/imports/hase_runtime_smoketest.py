#!/usr/bin/env python3

"""Verify that an installed HASEonGPU frontend owns its private runtime."""

import importlib.util
from pathlib import Path

import HASEonGPU
import pyInclude._runtime


module_path = Path(HASEonGPU.__file__).resolve()
runtime_path = Path(next(iter(pyInclude._runtime.__path__))).resolve()

print("module:", module_path)
print("runtime:", runtime_path)
assert "site-packages" in str(module_path) or "dist-packages" in str(module_path), module_path
assert (runtime_path / "calcPhiASE").is_file()
assert importlib.util.find_spec("HASEonGPU_Bindings") is None
