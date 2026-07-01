# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from _source_tree_import import ensure_hase_importable


scriptDir = Path(__file__).resolve().parent
defaultPhiAseConfigPath = scriptDir.parent / "config/hase-phiase.yaml"

ensure_hase_importable()

from HASEonGPU import (  # noqa: E402
    calcGainFromState,
    CrossSectionData,
    ExponentialEuler,
    GainMedium,
    MeshTopology,
    PhiASE,
    PumpProperties,
    Simulation,
    vtkWedge,
    writeParaviewState,
)


def printState(state):
    print(
        f"step={state.step:03d} "
        f"time={state.time:.3e}s "
        f"mean_beta={state.betaCells.mean():.6e} "
        f"mean_phi={state.phiAse.mean():.6e}"
    )


def writeVtkFields(state, vtkOutputDir=scriptDir, claddingAbsorption=1.0, crossSections=None, nTot=None):
    if state.phiAse is None:
        raise ValueError("VTK export requires state.phiAse")
    if crossSections is None:
        raise ValueError("VTK export requires crossSections for gain")
    if nTot is None:
        raise ValueError("VTK export requires nTot for gain")

    return vtkWedge(
        Path(vtkOutputDir) / f'laserPumpCladding_{state.step:03d}.vtk',
        state,
        fields={
            "betaCells": state.betaCells,
            "phiASE": state.phiAse,
            "dndtAse": state.dndtAse,
            "dndtPump": state.dndtPump,
            "cladAbs": state.phiAse * np.float64(claddingAbsorption),
            "localGain": calcGainFromState(state, crossSections, nTot),
        },
    )


def laserPumpCladdingMedium(numberOfLevels=10, thickness=None, cladAbsorption =5.5):
    materialPath = scriptDir / "data" / "pt.vtk"
    topology = MeshTopology.fromVtk(materialPath)

    if numberOfLevels is not None and topology.levels != numberOfLevels:
        raise ValueError(
            f"{materialPath} contains {topology.levels} levels, expected {numberOfLevels}"
        )
    if thickness is not None and not np.isclose(topology.thickness, thickness):
        raise ValueError(
            f"{materialPath} has thickness {topology.thickness}, expected {thickness}"
        )
    return GainMedium(topology=topology).withPhysicalProperties(
        betaCells=np.zeros((topology.numberOfPoints, topology.levels), dtype=np.float64),
        betaVolume=np.zeros((topology.numberOfTriangles, topology.levels - 1), dtype=np.float64),
        claddingCellTypes=np.zeros(topology.numberOfTriangles, dtype=np.uint32),
        refractiveIndices=np.asarray([1.83, 1.0, 1.83, 1.0], dtype=np.float32),
        reflectivities=np.zeros((topology.numberOfTriangles, 2), dtype=np.float32),
        nTot = 2.76e20,
        crystalTFluo = 9.5e-4,
        claddingNumber = 1,
        claddingAbsorption=cladAbsorption,
    )


def runExample(
    phiAseConfigPath=defaultPhiAseConfigPath,
    backend="UseConfig",
    timeSlices=150,
    # pumpSteps: pumped outer simulation steps; None pumps for all timeSlices.
    pumpSteps=50,
    vtkOutputDir=scriptDir,
    openPmdOutputDir=None,
    openpmdBackend="UseConfig",
    **AseOverride,
):
    vtkOutputDir = Path(vtkOutputDir)
    materialDir = scriptDir / "input"
    numberOfLevels = 10
    thickness = 0.7 / (numberOfLevels - 1)

    spectralProperties = CrossSectionData(
        wavelengthsAbsorption=np.loadtxt(materialDir / "lambda_a.txt"),
        crossSectionAbsorption=np.loadtxt(materialDir / "sigma_a.txt"),
        wavelengthsEmission=np.loadtxt(materialDir / "lambda_e.txt"),
        crossSectionEmission=np.loadtxt(materialDir / "sigma_e.txt"),
        resolution=np.loadtxt(materialDir / "lambda_a.txt").size,
    )

    pumpCrossSections = CrossSectionData.monochromatic(
        wavelength=940e-9,
        crossSectionAbsorption=0.778e-20,
        crossSectionEmission=0.195e-20,
    )
    absorption=5.5
    medium = laserPumpCladdingMedium(
        numberOfLevels=numberOfLevels,
        thickness=thickness,
        cladAbsorption=absorption
    )

    phiAse = PhiASE.fromYaml(
        phiAseConfigPath,
        spectralProperties=spectralProperties,
        **AseOverride
    )

    if backend != "UseConfig" : phiAse.backend=backend
    if openpmdBackend != "UseConfig" : phiAse.openpmdBackend=openpmdBackend


    pumpProperties=PumpProperties(
                         crossSections=pumpCrossSections,
                         intensity=16e3,
                         pumpDuration=1e-6,
                         pumpSubsteps=100,
                         temporaryFluorescence=1.0,
                         pumpSteps=pumpSteps,
                         wavelength=940e-9,
                         radiusX=1.5,
                         radiusY=1.5,
                         exponent=40,
                         backReflection=True,
                         reflectivity=1.0,
                         extraction=False)
    print(f"Running simulation with backend {phiAse.backend}")
    print(f"Using openPMD backend {phiAse.openpmdBackend}")
    simulation = Simulation(
        gainMedium=medium,
        pump=pumpProperties,
        phiASE=phiAse,
        timeIntegrationSolver=ExponentialEuler(),
        timeStep=2e-5,
        crossSections=spectralProperties,
    )
    simulation.onStep(printState)
    simulation.onStep(
        writeVtkFields,
        vtkOutputDir,
        absorption,
        spectralProperties,
        medium.get("nTot").value,
    )
    if openPmdOutputDir is not None:
        simulation.onStep(writeParaviewState, openPmdOutputDir, absorption)
    simulation.runSteps(timeSlices) # adjust this by number of steps
    return simulation.getLastState() # return the last state to confirm shape.

    # dndt_ASE, flux_clad
def main(argv=None):
    parser = argparse.ArgumentParser(description="Modern HASEonGPU laser-pump cladding example")
    parser.add_argument("--backend", type=str, default="UseConfig")
    parser.add_argument("--openpmd-backend", type=str, default="UseConfig")
    parser.add_argument("--timeSteps", type=int, default=150)
    parser.add_argument(
        "--pumpSteps",
        type=int,
        default=100,
        help=(
            "Number of outer simulation steps with pump contribution. "
            "Default: 100. Use a value matching --timeSteps to pump for the full run. "
            "This is distinct from "
            "PumpProperties.pumpSubsteps, which is the internal pump "
            "integration resolution."
        ),
    )
    parser.add_argument(
        "--phi-ase-config",
        type=Path,
        default=defaultPhiAseConfigPath,
        help="PhiASE run-control YAML. Defaults to config/hase-phiase.yaml.",
    )
    parser.add_argument("--vtk-output-dir", type=Path, default=scriptDir)
    parser.add_argument("--openpmd-output-dir", type=Path, default=None)
    parser.add_argument("--min-sample-range", type=int, default=None)
    parser.add_argument("--max-sample-range", type=int, default=None)
    parser.add_argument("--rng-seed", type=int, default=None)
    args = parser.parse_args(argv)

    aseOverrides = {}
    if args.min_sample_range is not None:
        aseOverrides["minSampleRange"] = args.min_sample_range
    if args.max_sample_range is not None:
        aseOverrides["maxSampleRange"] = args.max_sample_range
    if args.rng_seed is not None:
        aseOverrides["rngSeed"] = args.rng_seed

    state = runExample(
        args.phi_ase_config,
        args.backend,
        timeSlices=args.timeSteps,
        pumpSteps=args.pumpSteps,
        vtkOutputDir=args.vtk_output_dir,
        openPmdOutputDir=args.openpmd_output_dir,
        openpmdBackend=args.openpmd_backend,
        **aseOverrides,
    )
    print(f"phiAse shape: {state.phiAse.shape}")
    print(f"betaCells shape: {state.betaCells.shape}")


if __name__ == "__main__":
    main()
