# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Run laserPumpCladding with all model and run parameters supplied by YAML."""

from __future__ import annotations

import argparse
from pathlib import Path

from _source_tree_import import ensure_hase_importable


scriptDir = Path(__file__).resolve().parent
defaultConfigPath = scriptDir.parent / "config" / "laserPumpCladding.yaml"

ensure_hase_importable()

from HASEonGPU import Simulation  # noqa: E402
from laserPumpCladdingApi import printState, writeVtkFields  # noqa: E402
from pyInclude.openpmd.paraview import writeParaviewState  # noqa: E402


def buildSimulation(configPath=defaultConfigPath):
    """Construct the complete example without executing it."""
    return Simulation.from_yaml(configPath)


def runExample(
    configPath=defaultConfigPath,
    vtkOutputDir=scriptDir,
    openPmdOutputDir=None,
):
    """Run the YAML-defined simulation and return its final state."""
    simulation = buildSimulation(configPath)
    medium = simulation.gain_medium
    absorption = medium.get("claddingAbsorption").value

    print(f"Running simulation with backend {simulation.phi_ase.backend}")
    print(f"Using openPMD backend {simulation.phi_ase.openpmdBackend}")
    simulation.on_step(printState)
    simulation.on_step(
        writeVtkFields,
        Path(vtkOutputDir),
        absorption,
        simulation.cross_sections,
        medium.get("nTot").value,
    )
    if openPmdOutputDir is not None:
        simulation.on_step(writeParaviewState, Path(openPmdOutputDir), absorption)
    simulation.step()
    return simulation.get_last_state()


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "YAML-driven HASEonGPU laser-pump cladding example. Edit the YAML "
            "for model, solver, backend, and run parameters."
        )
    )
    parser.add_argument("--config", type=Path, default=defaultConfigPath)
    parser.add_argument("--vtk-output-dir", type=Path, default=scriptDir)
    parser.add_argument("--openpmd-output-dir", type=Path, default=None)
    args = parser.parse_args(argv)

    state = runExample(
        args.config,
        vtkOutputDir=args.vtk_output_dir,
        openPmdOutputDir=args.openpmd_output_dir,
    )
    print(f"phiAse shape: {state.phi_ase.shape}")
    print(f"betaVolume shape: {state.beta_volume.shape}")


if __name__ == "__main__":
    main()
