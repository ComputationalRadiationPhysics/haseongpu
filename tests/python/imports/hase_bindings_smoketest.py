#!/usr/bin/env python3

from pathlib import Path

import HASEonGPU
from HASEonGPU import ExperimentParameters, Result


module_path = Path(HASEonGPU.__file__).resolve()
module_path_str = str(module_path)

print("module:", module_path)
assert "site-packages" in module_path_str or "dist-packages" in module_path_str, module_path

result = Result()
params = ExperimentParameters()

print("Result type:", type(result))
print("phiAse entries:", result.num_phiAse)
print("ExperimentParameters type:", type(params))
print("minRaysPerSample:", params.minRaysPerSample)

assert type(result).__module__ == "HASEonGPU_Bindings.HASEonGPU"
assert type(params).__module__ == "HASEonGPU_Bindings.HASEonGPU"
assert params.minRaysPerSample == 100000
