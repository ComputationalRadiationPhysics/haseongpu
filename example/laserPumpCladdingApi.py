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

ensure_hase_importable()

from HASEonGPU import (  # noqa: E402
    calcGainFromState,
    CrossSectionData,
    FrozenPhiAseRungeKutta4,
    GainMedium,
    integrate_pump_profile,
    PlanarPumpRelay,
    PumpAngularDistribution,
    Pump,
    PumpSpectrum,
    SuperGaussianPumpProfile,
    SurfacePumpInjector,
    VolumeTopology,
    backendFlat,
    PhiASE,
    Simulation,
    SurfaceOptics,
    vtkWedge,
)
from pyInclude.openpmd.paraview import writeParaviewState  # noqa: E402


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
    backend="Host_Cpu_CpuOmpBlocks",
    spectralResolution=1000,
    **AseOverride,
):
    spectralProperties = laserPumpCladdingSpectralProperties(spectralResolution)
    medium = loadLaserPumpCladdingTet4Medium(materialPath)
    phiAse = PhiASE(
        spectralProperties=spectralProperties,
        propagationMode="forward",
        minRays=10000,
        maxRays=1000000,
        relativeStandardErrorThreshold=0.1,
        repetitions=4,
        adaptiveSteps=4,
        useReflections=True,
        reflectionMaxIterations=40,
        reflectionTolerance=1.0e-4,
        surfaceReservoirSize=32,
        monochromatic=False,
        backend=backend,
        openpmdBackend="auto",
        parallelMode="single",
        numDevices=1,
        nPerNode=1,
        **AseOverride,
    )
    phiAse.run(gainMedium=medium, crossSections=spectralProperties)
    return phiAse.getResults()


def printState(state):
    print(
        f"step={state.step:03d} "
        f"time={state.time:.3e}s "
        f"mean_beta={state.beta_volume.mean():.6e} "
        f"mean_phi={state.phi_ase.mean():.6e}"
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
    if state.phi_ase is None:
        raise ValueError("VTK export requires state.phi_ase")
    if crossSections is None:
        raise ValueError("VTK export requires crossSections for gain")
    if nTot is None:
        raise ValueError("VTK export requires nTot for gain")

    fields = {
        "betaVolume": state.beta_volume,
        "phiASE": state.phi_ase,
        "dndtAse": state.dndt_ase,
        "dndtPump": state.dndt_pump,
        "cladAbs": state.phi_ase * np.float64(claddingAbsorption),
        "localGain": calcGainFromState(state, crossSections, nTot),
    }
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
    return GainMedium(topology=topology).withPhysicalProperties(
        betaVolume=backendFlat(np.zeros(topology.numberOfCells, dtype=np.float64)),
        claddingCellTypes=np.zeros(topology.numberOfCells, dtype=np.uint32),
        nTot=2 * 1.388e20,
        crystalTFluo=9.41e-4,
        claddingNumber=1,
        claddingAbsorption=cladAbsorption,
    ).with_surface_optics(
        {
            "ase_bottom": SurfaceOptics(
                reflectivity=0.0,
                n_inside=1.83,
                n_outside=1.0,
            ),
            "ase_top": SurfaceOptics(
                reflectivity=0.0,
                n_inside=1.83,
                n_outside=1.0,
            ),
            "cladding": SurfaceOptics(
                reflectivity=0.0,
                n_inside=1.0,
                n_outside=1.0,
            ),
        }
    )


def buildSimulation(
    backend="Host_Cpu_CpuOmpBlocks",
    simulation_steps=150,
    pump_steps=50,
    ase_steps=150,
    openpmdBackend="auto",
    pre_pump=True,
    spectralResolution=1000,
    pumpRayCount=50000,
    pumpRngSeed=5489,
    reportTimings=False,
    outputSteps=None,
    **AseOverride,
):
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
    absorption = 5.5
    medium = laserPumpCladdingMedium(cladAbsorption=absorption)

    phiAseParameters = {
        "spectralProperties": spectralProperties,
        "propagationMode": "forward",
        "minRays": 10000,
        "maxRays": 1000000,
        "relativeStandardErrorThreshold": 0.1,
        "repetitions": 4,
        "adaptiveSteps": 4,
        "useReflections": True,
        "reflectionMaxIterations": 40,
        "reflectionTolerance": 1.0e-4,
        "surfaceReservoirSize": 32,
        "monochromatic": False,
        "backend": backend,
        "openpmdBackend": openpmdBackend,
        "parallelMode": "single",
        "numDevices": 1,
        "nPerNode": 1,
        "ase_steps": ase_steps,
    }
    phiAseParameters.update(AseOverride)
    phiAse = PhiASE(**phiAseParameters)


    pumpProfile = SuperGaussianPumpProfile(radius_u=1.5, radius_v=1.5, exponent=40)
    profileArea = (
        integrate_pump_profile(medium.topology, "ase_bottom", pumpProfile)
        if hasattr(medium.topology, "facePointIndices")
        else 1.0
    )
    pump = Pump(
        total_power=16e3 * profileArea,
        spectrum=PumpSpectrum.monochromatic(pumpWavelength),
        cross_sections=pumpCrossSections,
        ray_count=pumpRayCount,
        pump_steps=pump_steps,
        rng_seed=pumpRngSeed,
        angular_distribution=PumpAngularDistribution.collimated(),
        profile=pumpProfile,
    )
    return Simulation(
        gain_medium=medium,
        phi_ase=phiAse,
        time_integrator=FrozenPhiAseRungeKutta4(),
        time_step_size=2e-5,
        simulation_steps=simulation_steps,
        cross_sections=spectralProperties,
        pre_pump=pre_pump,
        report_timings=reportTimings,
        output_steps=None if outputSteps is None else tuple(int(step) for step in outputSteps),
    ).add_pump(
        pump,
        injection_method=SurfacePumpInjector(surface_domains="ase_bottom"),
        relays=(PlanarPumpRelay.retroreflect("ase_top"),),
    )


def runExample(
    backend="Host_Cpu_CpuOmpBlocks",
    simulation_steps=150,
    pump_steps=50,
    ase_steps=150,
    vtkOutputDir=scriptDir,
    openPmdOutputDir=None,
    openpmdBackend="auto",
    pre_pump=True,
    spectralResolution=1000,
    pumpRayCount=50000,
    pumpRngSeed=5489,
    reportTimings=False,
    outputSteps=None,
    **AseOverride,
):
    simulation = buildSimulation(
        backend=backend,
        simulation_steps=simulation_steps,
        pump_steps=pump_steps,
        ase_steps=ase_steps,
        openpmdBackend=openpmdBackend,
        pre_pump=pre_pump,
        spectralResolution=spectralResolution,
        pumpRayCount=pumpRayCount,
        pumpRngSeed=pumpRngSeed,
        reportTimings=reportTimings,
        outputSteps=outputSteps,
        **AseOverride,
    )
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
        simulation.on_step(writeParaviewState, openPmdOutputDir, absorption)
    simulation.step()
    return simulation.get_last_state()  # return the last state to confirm shape.

    # dndt_ASE, flux_clad
def main(argv=None):
    parser = argparse.ArgumentParser(description="Inline-API HASEonGPU laser-pump cladding example")
    parser.add_argument("--backend", type=str, default="Host_Cpu_CpuOmpBlocks")
    parser.add_argument("--openpmd-backend", type=str, default="auto")
    parser.add_argument("--simulation-steps", type=int, default=150)
    parser.add_argument(
        "--output-steps",
        type=int,
        nargs="+",
        default=None,
        help=(
            "Completed one-based step indices to emit. By default, emit every "
            "completed step."
        ),
    )
    parser.add_argument(
        "--pump-steps",
        type=int,
        default=50,
        help=(
            "Number of outer simulation steps with pump contribution. "
            "Default: 50. Zero disables this pump. Use a value matching "
            "--simulation-steps to pump for the full run."
        ),
    )
    parser.add_argument(
        "--ase-steps",
        type=int,
        default=150,
        help="Number of outer steps with ASE. Zero disables ASE. Default: 150.",
    )
    parser.add_argument("--vtk-output-dir", type=Path, default=scriptDir)
    parser.add_argument("--openpmd-output-dir", type=Path, default=None)
    parser.add_argument(
        "--disable-pre-pump",
        action="store_true",
        help="Run ASE during the first pump time step instead of seeding beta without ASE.",
    )
    parser.add_argument("--tet4-input", type=Path, default=None)
    parser.add_argument("--phiase-only", action="store_true")
    parser.add_argument("--rng-seed", type=int, default=None)
    parser.add_argument(
        "--pump-ray-count", type=int, default=50000,
        help="Equal-power launch rays per pump source. Default: 50000.",
    )
    parser.add_argument("--pump-rng-seed", type=int, default=5489)
    parser.add_argument(
        "--spectral-resolution",
        type=int,
        default=1000,
        help="Backend spectral interpolation resolution. Default: 1000.",
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
    if args.phiase_only:
        if args.tet4_input is None:
            parser.error("--phiase-only requires --tet4-input")
        aseOverrides.setdefault("propagationMode", "forward")
        result = runTet4PhiAseInput(
            args.tet4_input,
            args.backend,
            args.spectral_resolution,
            **aseOverrides,
        )
        phi = np.asarray(result.phiAse)
        print(f"phiAse shape: {phi.shape}")
        print(f"meanPhi: {float(phi.mean()):.17g}")
        return

    state = runExample(
        args.backend,
        simulation_steps=args.simulation_steps,
        pump_steps=args.pump_steps,
        ase_steps=args.ase_steps,
        vtkOutputDir=args.vtk_output_dir,
        openPmdOutputDir=args.openpmd_output_dir,
        openpmdBackend=args.openpmd_backend,
        pre_pump=not args.disable_pre_pump,
        spectralResolution=args.spectral_resolution,
        reportTimings=args.timings,
        pumpRayCount=args.pump_ray_count,
        pumpRngSeed=args.pump_rng_seed,
        outputSteps=args.output_steps,
        **aseOverrides,
    )
    print(f"phiAse shape: {state.phi_ase.shape}")
    print(f"betaVolume shape: {state.beta_volume.shape}")


if __name__ == "__main__":
    main()
