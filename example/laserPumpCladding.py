# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np

from _source_tree_import import ensure_hase_importable


scriptDir = Path(__file__).resolve().parent
defaultPhiAseConfigPath = Path(
    os.environ.get(
        "HASE_PHIASE_CONFIG",
        scriptDir.parent / "config/hase-phiase.yaml",
    )
)

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
    SurfaceOptics,
    vtkWedge,
    writeParaviewState,
)


def _loadLaserPumpCladdingRawSpectra():
    materialDir = scriptDir / "input"
    return (
        np.loadtxt(materialDir / "lambda_a.txt"),
        np.loadtxt(materialDir / "sigma_a.txt"),
        np.loadtxt(materialDir / "lambda_e.txt"),
        np.loadtxt(materialDir / "sigma_e.txt"),
    )


def laserPumpCladdingSpectralProperties(spectralResolution=1000):
    """Return the raw material spectrum; transport resampling belongs to the backend."""
    (
        raw_wavelengths_absorption,
        raw_cross_section_absorption,
        raw_wavelengths_emission,
        raw_cross_section_emission,
    ) = _loadLaserPumpCladdingRawSpectra()
    return CrossSectionData(
        wavelengthsAbsorption=raw_wavelengths_absorption,
        crossSectionAbsorption=raw_cross_section_absorption,
        wavelengthsEmission=raw_wavelengths_emission,
        crossSectionEmission=raw_cross_section_emission,
        resolution=spectralResolution,
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
    spectralResolution=1000,
    **AseOverride,
):
    spectralProperties = laserPumpCladdingSpectralProperties(spectralResolution)
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
    if getattr(state, "volumePhiAse", None) is not None:
        fields["volumePhiASE"] = state.volumePhiAse
    if getattr(state, "volumeDndtAse", None) is not None:
        fields["volumeDndtAse"] = state.volumeDndtAse
    path = Path(vtkOutputDir) / f"laserPumpCladding_{state.step:03d}.vtk"
    if hasattr(state.topology, "cellPointIndices"):
        return _writeTet4StateVtk(path, state, fields)
    return vtkWedge(path, state, fields=fields)


BOTTOM_ASE_SURFACE_ID = 1
TOP_ASE_SURFACE_ID = 2
CLADDING_SURFACE_ID = 3
NUMBER_OF_Z_LAYERS = 10


def _assignLegacyTet4SurfaceDomains(topology):
    """Attach the legacy optical regions to its geometrically identical Tet4 mesh."""
    sample_points = np.asarray(topology.samplePoints, dtype=np.float64).copy()
    points = np.asarray(topology.points, dtype=np.float64)
    exterior = topology.neighborCells < 0
    z = points[:, 2]
    face_z = z[np.asarray(topology.facePointIndices, dtype=np.uint32)]
    bottom = exterior & np.all(np.isclose(face_z, np.min(z)), axis=2)
    top = exterior & np.all(np.isclose(face_z, np.max(z)), axis=2)
    side = exterior & ~bottom & ~top
    if not (np.any(bottom) and np.any(top) and np.any(side)):
        raise ValueError(
            "ptTet4.vtk must contain bottom, top, and cladding exterior faces"
        )

    topology = topology.withCellDomains(
        domain=1,
        name="gain_medium",
        where="all",
    ).withSurfaceDomains(
        [
            {
                "domain": BOTTOM_ASE_SURFACE_ID,
                "name": "ase_bottom",
                "faceIndices": np.argwhere(bottom),
            },
            {
                "domain": TOP_ASE_SURFACE_ID,
                "name": "ase_top",
                "faceIndices": np.argwhere(top),
            },
            {
                "domain": CLADDING_SURFACE_ID,
                "name": "cladding",
                "faceIndices": np.argwhere(side),
            },
        ]
    )
    topology.samplePoints = sample_points
    return topology


def laserPumpCladdingMedium(cladAbsorption=5.5):
    materialPath = scriptDir / "data" / "ptTet4.vtk"
    topology = _assignLegacyTet4SurfaceDomains(VolumeTopology.fromVtk(materialPath))
    refractiveIndices = np.asarray([1.83, 1.0, 1.83, 1.0], dtype=np.float32)
    return GainMedium(topology=topology).withPhysicalProperties(
        betaCells=backendFlat(np.zeros(topology.numberOfSamplePoints, dtype=np.float64)),
        betaVolume=backendFlat(np.zeros(topology.numberOfCells, dtype=np.float64)),
        claddingCellTypes=np.zeros(topology.numberOfCells, dtype=np.uint32),
        refractiveIndices=refractiveIndices,
        reflectivities=np.zeros((topology.numberOfCells, 2), dtype=np.float32),
        nTot=2 * 1.388e20,
        crystalTFluo=9.41e-4,
        claddingNumber=1,
        claddingAbsorption=cladAbsorption,
    ).withSurfaceOptics(
        {
            "ase_bottom": SurfaceOptics(
                reflectivity=0.0,
                n_inside=refractiveIndices[0],
                n_outside=refractiveIndices[1],
            ),
            "ase_top": SurfaceOptics(
                reflectivity=0.0,
                n_inside=refractiveIndices[2],
                n_outside=refractiveIndices[3],
            ),
            "cladding": SurfaceOptics(
                reflectivity=0.0,
                n_inside=refractiveIndices[1],
                n_outside=refractiveIndices[1],
            ),
        }
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
    prePump=True,
    spectralResolution=1000,
    reportTimings=False,
    **AseOverride,
):
    vtkOutputDir = Path(vtkOutputDir)
    spectralProperties = laserPumpCladdingSpectralProperties(spectralResolution)

    pumpWavelength = 940e-9
    (
        raw_wavelengths_absorption,
        raw_cross_section_absorption,
        raw_wavelengths_emission,
        raw_cross_section_emission,
    ) = _loadLaserPumpCladdingRawSpectra()
    pumpCrossSections = CrossSectionData.monochromatic(
        wavelength=pumpWavelength,
        crossSectionAbsorption=np.interp(
            pumpWavelength * 1.0e9,
            raw_wavelengths_absorption,
            raw_cross_section_absorption,
        ),
        crossSectionEmission=np.interp(
            pumpWavelength * 1.0e9,
            raw_wavelengths_emission,
            raw_cross_section_emission,
        ),
    )
    absorption=5.5
    medium = laserPumpCladdingMedium(cladAbsorption=absorption)

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
                         wavelength=pumpWavelength,
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
        prePump=prePump,
        reportTimings=reportTimings,
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
    parser.add_argument(
        "--disable-pre-pump",
        action="store_true",
        help="Run ASE during the first pump time step instead of seeding beta without ASE.",
    )
    parser.add_argument("--tet4-input", type=Path, default=None)
    parser.add_argument("--phiase-only", action="store_true")
    parser.add_argument("--rng-seed", type=int, default=None)
    parser.add_argument(
        "--spectral-resolution",
        type=int,
        default=1000,
        help="Backend spectral interpolation resolution. Default: 1000.",
    )
    parser.add_argument(
        "--disable-reflections",
        action="store_true",
        help="Disable ASE surface reflections.",
    )
    parser.add_argument(
        "--timings",
        action="store_true",
        help="Print frontend timing split for compiled transport, snapshots, and callbacks.",
    )
    args = parser.parse_args(argv)

    aseOverrides = {}
    if args.rng_seed is not None:
        aseOverrides["rngSeed"] = args.rng_seed
    if args.disable_reflections:
        aseOverrides["useReflections"] = False

    if args.phiase_only:
        if args.tet4_input is None:
            parser.error("--phiase-only requires --tet4-input")
        aseOverrides.setdefault("propagationMode", "forward")
        result = runTet4PhiAseInput(
            args.tet4_input,
            args.phi_ase_config,
            args.backend,
            args.spectral_resolution,
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
        prePump=not args.disable_pre_pump,
        spectralResolution=args.spectral_resolution,
        reportTimings=args.timings,
        **aseOverrides,
    )
    print(f"phiAse shape: {state.phiAse.shape}")
    print(f"betaCells shape: {state.betaCells.shape}")


if __name__ == "__main__":
    main()
