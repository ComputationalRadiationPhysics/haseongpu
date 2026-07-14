#!/usr/bin/env python3
# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

from pathlib import Path

import HASEonGPU
from HASEonGPU import GainMedium, Grid, MeshTopology, PhiASE, TransportResult


module_path = Path(HASEonGPU.__file__).resolve()
module_path_str = str(module_path)

print("module:", module_path)
assert "site-packages" in module_path_str or "dist-packages" in module_path_str, module_path

topology = MeshTopology.fromGrid(Grid(xExtent=1, yExtent=1, zExtent=1, tileSizeZ=0.5))
gain_medium = GainMedium(topology=topology)
phi_ase = PhiASE(backend="Host_Cpu_CpuSerial")
result = TransportResult()

print("GainMedium type:", type(gain_medium))
print("PhiASE type:", type(phi_ase))
print("TransportResult type:", type(result))
print("numberOfPrisms:", gain_medium.numberOfPrisms)

assert gain_medium.numberOfPrisms == 4
assert result.phiAse == []
assert phi_ase.minRays == 100000

for legacy_name in ("HostMesh", "ExperimentParameters", "ComputeParameters", "Mesh"):
    assert not hasattr(HASEonGPU, legacy_name), legacy_name
