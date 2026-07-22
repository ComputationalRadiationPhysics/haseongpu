#!/usr/bin/env python3

"""Verify that an installed thin frontend resolves its durable runtime."""

import importlib.util
from pathlib import Path

import HASEonGPU
import openpmd_api
import pyInclude._runtime
from pyInclude._runtime import runtime_config, runtime_root


module_path = Path(HASEonGPU.__file__).resolve()
runtime_path = Path(next(iter(pyInclude._runtime.__path__))).resolve()
native_runtime = runtime_root()
provider_path = str(runtime_config().HASE_OPENPMD_PYTHON_PACKAGE_DIR or "")
openpmd_path = Path(openpmd_api.__file__).resolve()

print("module:", module_path)
print("frontend metadata:", runtime_path)
print("native runtime:", native_runtime)
print("openPMD Python provider:", openpmd_path)
assert "site-packages" in str(module_path) or "dist-packages" in str(module_path), module_path
assert not (runtime_path / "calcPhiASE").exists()
assert (native_runtime / "calcPhiASE").is_file() or (
    native_runtime / "python" / "pyInclude" / "_runtime" / "calcPhiASE"
).is_file()
assert Path(runtime_config().HASE_RUNTIME_DIR).resolve() == native_runtime
if provider_path:
    openpmd_path.relative_to(Path(provider_path).resolve())
assert importlib.util.find_spec("HASEonGPU_Bindings") is None
