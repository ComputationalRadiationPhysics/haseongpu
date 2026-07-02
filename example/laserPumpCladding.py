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
    FrozenPhiAseRungeKutta4,
    GainMedium,
    VolumeTopology,
    backendFlat,
    PhiASE,
    PumpProperties,
    Simulation,
    vtkWedge,
    writeParaviewState,
)


def laserPumpCladdingSpectralProperties():
    materialDir = scriptDir / "input"
    return CrossSectionData(
        wavelengthsAbsorption=np.loadtxt(materialDir / "lambda_a.txt"),
        crossSectionAbsorption=np.loadtxt(materialDir / "sigma_a.txt"),
        wavelengthsEmission=np.loadtxt(materialDir / "lambda_e.txt"),
        crossSectionEmission=np.loadtxt(materialDir / "sigma_e.txt"),
        resolution=np.loadtxt(materialDir / "lambda_a.txt").size,
    )


def loadLaserPumpCladdingTet4Medium(materialPath):
    """Load a Tet4 laserPumpCladding state for PhiASE-only runs.

    Converted legacy VTK fixtures store both point and cell pump data.  The
    forward openPMD backend writes PhiASE results on tetrahedral cells.
    """
    return GainMedium.fromVtk(materialPath)


def runTet4PhiAseInput(
    materialPath,
    phiAseConfigPath=defaultPhiAseConfigPath,
    backend="UseConfig",
    **AseOverride,
):
    spectralProperties = laserPumpCladdingSpectralProperties()
    medium = loadLaserPumpCladdingTet4Medium(materialPath)
    phiAse = PhiASE.fromYaml(
        phiAseConfigPath,
        spectralProperties=spectralProperties,
        **AseOverride,
    )
    if backend != "UseConfig":
        phiAse.backend = backend
    phiAse.run(gainMedium=medium, crossSections=spectralProperties)
    return phiAse.getResults()


def printState(state):
    print(
        f"step={state.step:03d} "
        f"time={state.time:.3e}s "
        f"mean_beta={state.betaCells.mean():.6e} "
        f"mean_phi={state.phiAse.mean():.6e}"
    )


def _writeScalarArray(handle, name, values, count):
    arr = np.asarray(values).reshape(-1, order="F")
    if arr.size != count:
        raise ValueError(f"{name} has {arr.size} values, expected {count}")
    handle.write(f"SCALARS {name} double 1\n")
    handle.write("LOOKUP_TABLE default\n")
    for value in arr:
        handle.write(f"{float(value):.17g}\n")


def _writeTet4StateVtk(path, state, fields):
    topology = state.topology
    points = np.asarray(topology.points, dtype=np.float64)
    cells = np.asarray(topology.cellPointIndices, dtype=np.uint32)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    point_count = points.shape[0]
    cell_count = cells.shape[0]
    point_fields = {name: value for name, value in fields.items() if np.asarray(value).size == point_count}
    cell_fields = {name: value for name, value in fields.items() if np.asarray(value).size == cell_count}
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# vtk DataFile Version 2.0\n")
        handle.write("HASEonGPU laserPumpCladding Tet4 state\n")
        handle.write("ASCII\n")
        handle.write("DATASET UNSTRUCTURED_GRID\n")
        handle.write(f"POINTS {point_count} double\n")
        for x, y, z in points:
            handle.write(f"{x:.17g} {y:.17g} {z:.17g}\n")
        handle.write(f"CELLS {cell_count} {cell_count * 5}\n")
        for cell in cells:
            handle.write("4 " + " ".join(str(int(vertex)) for vertex in cell) + "\n")
        handle.write(f"CELL_TYPES {cell_count}\n")
        handle.write(("10\n" * cell_count))
        if point_fields:
            handle.write(f"POINT_DATA {point_count}\n")
            for name, values in point_fields.items():
                _writeScalarArray(handle, name, values, point_count)
        if cell_fields:
            handle.write(f"CELL_DATA {cell_count}\n")
            for name, values in cell_fields.items():
                _writeScalarArray(handle, name, values, cell_count)
    return path


def writeVtkFields(state, vtkOutputDir=scriptDir, claddingAbsorption=1.0, crossSections=None, nTot=None):
    if state.phiAse is None:
        raise ValueError("VTK export requires state.phiAse")
    if crossSections is None:
        raise ValueError("VTK export requires crossSections for gain")
    if nTot is None:
        raise ValueError("VTK export requires nTot for gain")

    fields = {
        "betaCells": state.betaCells,
        "betaVolume": state.betaVolume,
        "phiASE": state.phiAse,
        "dndtAse": state.dndtAse,
        "dndtPump": state.dndtPump,
        "cladAbs": state.phiAse * np.float64(claddingAbsorption),
        "localGain": calcGainFromState(state, crossSections, nTot),
    }
    path = Path(vtkOutputDir) / f"laserPumpCladding_{state.step:03d}.vtk"
    if hasattr(state.topology, "cellPointIndices"):
        return _writeTet4StateVtk(path, state, fields)
    return vtkWedge(path, state, fields=fields)



def laserPumpCladdingMedium(numberOfLevels=10, thickness=None, cladAbsorption =5.5):
    materialPath = scriptDir / "data" / "pt.vtk"
    topology = VolumeTopology.fromVtk(materialPath)

    if numberOfLevels is not None and topology.structuredNumberOfLevels != numberOfLevels:
        raise ValueError(
            f"{materialPath} contains {topology.structuredNumberOfLevels} levels, expected {numberOfLevels}"
        )
    if thickness is not None and not np.isclose(topology.structuredThickness, thickness):
        raise ValueError(
            f"{materialPath} has thickness {topology.structuredThickness}, expected {thickness}"
        )
    return GainMedium(topology=topology).withPhysicalProperties(
        betaCells=backendFlat(np.zeros(topology.numberOfSamplePoints, dtype=np.float64)),
        betaVolume=backendFlat(np.zeros(topology.numberOfCells, dtype=np.float64)),
        claddingCellTypes=np.zeros(topology.numberOfCells, dtype=np.uint32),
        refractiveIndices=np.asarray([1.83, 1.0, 1.83, 1.0], dtype=np.float32),
        reflectivities=np.zeros((topology.numberOfCells, 2), dtype=np.float32),
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
    openPmdOutputDir=None,
    openpmdBackend="UseConfig",
    enableASE=True,
    **AseOverride,
):
    vtkOutputDir = Path(vtkOutputDir)
    numberOfLevels = 10
    thickness = 0.7 / (numberOfLevels - 1)

    spectralProperties = laserPumpCladdingSpectralProperties()

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

    AseOverride.setdefault("forwardRayLength", 1.0)
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
        timeIntegrationSolver=FrozenPhiAseRungeKutta4(),
        timeStep=2e-5,
        crossSections=spectralProperties,
        enableASE=enableASE,
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
    parser.add_argument(
        "--disable-ase",
        action="store_true",
        help="Disable ASE depletion during the time-stepped pump simulation.",
    )
    parser.add_argument("--tet4-input", type=Path, default=None)
    parser.add_argument("--phiase-only", action="store_true")
    parser.add_argument("--min-sample-range", type=int, default=None)
    parser.add_argument("--max-sample-range", type=int, default=None)
    parser.add_argument("--rng-seed", type=int, default=None)
    parser.add_argument(
        "--disable-reflections",
        action="store_true",
        help="Disable ASE surface reflections; required by the current forward Tet4 backend.",
    )
    args = parser.parse_args(argv)

    aseOverrides = {}
    if args.min_sample_range is not None:
        aseOverrides["minSampleRange"] = args.min_sample_range
    if args.max_sample_range is not None:
        aseOverrides["maxSampleRange"] = args.max_sample_range
    if args.rng_seed is not None:
        aseOverrides["rngSeed"] = args.rng_seed
    if args.disable_reflections:
        aseOverrides["useReflections"] = False

    if args.phiase_only:
        if args.tet4_input is None:
            parser.error("--phiase-only requires --tet4-input")
        aseOverrides.setdefault("propagationMode", "forward")
        aseOverrides.setdefault("forwardRayLength", 1.0)
        result = runTet4PhiAseInput(
            args.tet4_input,
            args.phi_ase_config,
            args.backend,
            **aseOverrides,
        )
        phi = np.asarray(result.phiAse)
        print(f"phiAse shape: {phi.shape}")
        print(f"meanPhi: {float(phi.mean()):.17g}")
        return

    state = runExample(
        args.phi_ase_config,
        args.backend,
        timeSlices=args.timeSteps,
        pumpSteps=args.pumpSteps,
        vtkOutputDir=args.vtk_output_dir,
        openPmdOutputDir=args.openpmd_output_dir,
        openpmdBackend=args.openpmd_backend,
        enableASE=not args.disable_ase,
        **aseOverrides,
    )
    print(f"phiAse shape: {state.phiAse.shape}")
    print(f"betaCells shape: {state.betaCells.shape}")


if __name__ == "__main__":
    main()
