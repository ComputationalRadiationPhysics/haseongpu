# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from pyInclude import FrozenPhiAseRungeKutta4

scriptDir = Path(__file__).resolve().parent
repoRoot = scriptDir.parents[1]
buildPythonRoot = repoRoot / "build" / "python"
defaultPhiAseConfigPath = scriptDir / "config/phiASE.yaml"

sys.path.insert(0, str(repoRoot))
if buildPythonRoot.is_dir():
    sys.path.insert(0, str(buildPythonRoot))

from HASEonGPU import (  # noqa: E402
    OneDimensionalZTraversal,
    calcGainFromState,
    Constants,
    CrossSectionData,
    GainMedium,
    MeshTopology,
    PhiASE,
    PumpProperties,
    Simulation,
    PumpRadiationProfile,
    RungeKutta4,
    vtkWedge,
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
    from scipy.io import loadmat

    materialPath = scriptDir / "legacy" / "pt.mat"
    thickness = 0.7 / (numberOfLevels - 1) if thickness is None else thickness
    material = loadmat(materialPath)
    points = np.asarray(material["p"], dtype=np.float64)
    triangles = np.asarray(material["t"], dtype=np.int64) - 1
    if np.any(triangles < 0):
        raise ValueError(f"{materialPath} contains non-MATLAB-style triangle indices")

    topology = MeshTopology(
        points=points[:, :2],
        trianglePointIndices=triangles.astype(np.uint32),
        levels=numberOfLevels,
        thickness=thickness,
        metadata={"source": str(materialPath), "format": "matlab-pt"},
    )
    return GainMedium(topology=topology).withPhysicalProperties(
        betaCells=np.zeros((topology.numberOfPoints, topology.levels), dtype=np.float64),
        betaVolume=np.zeros((topology.numberOfTriangles, topology.levels - 1), dtype=np.float64),
        claddingCellTypes=np.zeros(topology.numberOfTriangles, dtype=np.uint32),
        refractiveIndices=np.asarray([1.83, 1.0, 1.83, 1.0], dtype=np.float32),
        reflectivities=np.zeros((2, topology.numberOfTriangles), dtype=np.float32),
        nTot=2 * 1.388e20,
        crystalTFluo=9.41e-4,
        claddingNumber=1,
        claddingAbsorption=cladAbsorption,
    )


def runExample(
    phiAseConfigPath=defaultPhiAseConfigPath,
    backend="UseConfig",
    timeSlices=150,
    # pumpSteps: pumped outer simulation steps; None pumps for all timeSlices.
    pumpSteps=50,
    vtkOutputDir=scriptDir,
    enableAse=True,
    prePump=True,
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
        resolution=1000,
    )

    pumpProfile = PumpRadiationProfile(
        intensity=16e3, # KW/cm^2
        wavelengths=[940e-9],
        waist=(1.5, 1.5),
        superGaussianOrder=40,
        propagationDirection=(0.0, 0.0, 1.0),
        backReflection=True,
        reflectivity=1.0,
    )
    absorption=5.5 # constant absoption coefficient
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


    pumpProperties=PumpProperties(
                         crossSections=spectralProperties,
                         profile=pumpProfile,
                         pumpSteps=pumpSteps, # allows stopping pumping earlier to see gain degrading
                         solver=OneDimensionalZTraversal(),
                         extraction=False)
    print(f"Running simulation with backend {phiAse.backend}")
    simulation = Simulation(
        gainMedium=medium,
        pump=pumpProperties,
        phiASE=phiAse,
        timeIntegrationSolver=FrozenPhiAseRungeKutta4(),
        timeStep=2e-5, #20 μs
        crossSections=spectralProperties,
        enableAse=enableAse,
        prePump=prePump,
    )

    simulation.onStep(printState)
    simulation.onStep(
        writeVtkFields,
        vtkOutputDir,
        absorption,
        spectralProperties,
        medium.get("nTot").value,
    )
    simulation.runSteps(timeSlices) # adjust this by number of steps
    return simulation.getLastState() # return the last state to confirm shape.

    # dndt_ASE, flux_clad
def main(argv=None):
    parser = argparse.ArgumentParser(description="Modern HASEonGPU laser-pump cladding example")
    parser.add_argument("--backend", type=str, default="UseConfig")
    parser.add_argument("--timeSteps", type=int, default=150)
    parser.add_argument(
        "--pumpSteps",
        type=int,
        default=100,
        help=(
            "Number of outer simulation steps with pump contribution. "
            "Default: 100. Use a value matching --timeSteps to pump for the full run. "
            "The continuous pump solver evaluates dndtPump once per outer step; "
            "there are no internal pump substeps."
        ),
    )
    parser.add_argument("--phi-ase-config", type=Path, default=defaultPhiAseConfigPath)
    parser.add_argument("--vtk-output-dir", type=Path, default=scriptDir)
    parser.add_argument("--disable-ase", action="store_true", help="Skip PhiASE backend runs and use zero ASE depletion")
    parser.add_argument(
        "--disable-pre-pump",
        action="store_true",
        help="Run ASE from the first simulation step instead of first seeding beta with pump only",
    )
    args = parser.parse_args(argv)

    state = runExample(
        args.phi_ase_config,
        args.backend,
        timeSlices=args.timeSteps,
        pumpSteps=args.pumpSteps,
        vtkOutputDir=args.vtk_output_dir,
        enableAse=not args.disable_ase,
        prePump=not args.disable_pre_pump,
    )
    print(f"phiAse shape: {state.phiAse.shape}")
    print(f"betaCells shape: {state.betaCells.shape}")


if __name__ == "__main__":
    main()
