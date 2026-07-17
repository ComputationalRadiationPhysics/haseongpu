# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

from pathlib import Path
import platform
import sys

import HASEonGPU
import HASEonGPU_Bindings
import HASEonGPU_Bindings.HASEonGPU as native_bindings

sys.path.insert(0, "/src/example/python_example")
import laserPumpCladding


print(f"python={platform.python_version()}", flush=True)
print(f"HASEonGPU={Path(HASEonGPU.__file__).resolve()}", flush=True)
print(f"HASEonGPU_Bindings={Path(HASEonGPU_Bindings.__file__).resolve()}", flush=True)
print(f"native_bindings={Path(native_bindings.__file__).resolve()}", flush=True)

state = laserPumpCladding.runExample(
    phiAseConfigPath=Path("/generation/phiASE-no-reflections.yaml"),
    backend="Host_Cpu_CpuOmpBlocks",
    timeSlices=6,
    pumpSteps=3,
    vtkOutputDir=Path("/output"),
    enableAse=True,
    prePump=True,
    rngSeed=5489,
)

print(f"final_step={state.step}", flush=True)
print(f"final_phi_shape={state.phiAse.shape}", flush=True)
