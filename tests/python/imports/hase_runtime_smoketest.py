#!/usr/bin/env python3

"""Verify that an installed thin frontend resolves its durable runtime."""

import importlib.util
from pathlib import Path

import HASEonGPU
import pyInclude._runtime
from pyInclude._runtime import runtime_config, runtime_root


module_path = Path(HASEonGPU.__file__).resolve()
runtime_path = Path(next(iter(pyInclude._runtime.__path__))).resolve()
native_runtime = runtime_root()

print("module:", module_path)
print("frontend metadata:", runtime_path)
print("native runtime:", native_runtime)
assert "site-packages" in str(module_path) or "dist-packages" in str(module_path), module_path
assert not (runtime_path / "calcPhiASE").exists()
assert (native_runtime / "calcPhiASE").is_file() or (
    native_runtime / "python" / "pyInclude" / "_runtime" / "calcPhiASE"
).is_file()
assert Path(runtime_config().HASE_RUNTIME_DIR).resolve() == native_runtime
assert importlib.util.find_spec("HASEonGPU_Bindings") is None
