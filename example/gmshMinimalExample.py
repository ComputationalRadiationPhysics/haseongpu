# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

from pathlib import Path
from tempfile import TemporaryDirectory

import gmsh
import numpy as np

from _source_tree_import import ensure_hase_importable

ensure_hase_importable()

from HASEonGPU import (
    GainMedium,
    MeshTopology,
    PhiASE,
    PumpProperties,
    RungeKutta4,
    Simulation,
    SpectralDecomposition,
    vtkWedge,
)


def initFunc(simulation):
    medium = simulation.gainMedium
    print(f"gmsh topology: {medium.numberOfTriangles} triangles, {medium.numberOfPrisms} prisms")


def printState(state):
    print(
        f"step={state.step:03d} "
        f"time={state.time:.3e}s "
        f"mean_beta={state.betaCells.mean():.6e} "
        f"mean_phi={state.phiAse.mean():.6e}"
    )


def writeVtkState(state, outputFile):
    vtkWedge(outputFile, state)


def cylindrical_core_cladding_surfaces(core_radius, cladding_radius, *, mesh_size):
    geo = gmsh.model.geo
    center = geo.addPoint(0.0, 0.0, 0.0, mesh_size)
    rings = []
    for radius in (core_radius, cladding_radius):
        rings.append(
            [
                geo.addPoint(radius, 0.0, 0.0, mesh_size),
                geo.addPoint(0.0, radius, 0.0, mesh_size),
                geo.addPoint(-radius, 0.0, 0.0, mesh_size),
                geo.addPoint(0.0, -radius, 0.0, mesh_size),
            ]
        )
    core_points, cladding_points = rings
    core_arcs = [geo.addCircleArc(core_points[i], center, core_points[(i + 1) % 4]) for i in range(4)]
    cladding_arcs = [geo.addCircleArc(cladding_points[i], center, cladding_points[(i + 1) % 4]) for i in range(4)]
    core = geo.addPlaneSurface([geo.addCurveLoop(core_arcs)])
    cladding = geo.addPlaneSurface([geo.addCurveLoop(cladding_arcs), geo.addCurveLoop([-arc for arc in core_arcs])])
    return core, cladding


def write_minimal_gmsh_mesh(filename):
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("hase_cylindrical_core_cladding")
        core, cladding = cylindrical_core_cladding_surfaces(0.5, 0.8, mesh_size=0.4)
        gmsh.model.geo.synchronize()

        gmsh.model.addPhysicalGroup(2, [core], 20)
        gmsh.model.setPhysicalName(2, 20, "Core")
        gmsh.model.addPhysicalGroup(2, [cladding], 21)
        gmsh.model.setPhysicalName(2, 21, "CladdingShell")

        gmsh.model.mesh.generate(2)
        gmsh.write(str(filename))
    finally:
        gmsh.finalize()


def main():
    with TemporaryDirectory() as tmpdir:
        gmsh_file = Path(tmpdir) / "minimal_core_cladding.msh"
        write_minimal_gmsh_mesh(gmsh_file)

        topology = MeshTopology.fromFile(gmsh_file, format="gmsh", numberOfLevels=6, thickness=0.25)
        medium = GainMedium(topology=topology)

        print(medium.get("betaCells").expectedShape)
        print("gmsh claddingCellTypes:", medium.get("claddingCellTypes").value)

        medium.get("betaCells").value = np.zeros(medium.get("betaCells").expectedShape)
        medium.get("betaVolume").value = np.zeros(medium.get("betaVolume").expectedShape)
        medium.get("claddingCellTypes").value = np.asarray(
            medium.get("claddingCellTypes").value, dtype=np.uint32
        ).reshape(medium.get("claddingCellTypes").expectedShape)
        medium.get("refractiveIndices").value = np.asarray([2.0, 1.0, 3.0, 4.0], dtype=np.float32)
        medium.get("reflectivities").value = np.zeros(
            medium.get("reflectivities").expectedShape, dtype=np.float32
        )
        medium.get("nTot").value = 1.388e20 * 2.0  # Doping density [1/cm^3]
        medium.get("crystalTFluo").value = 9.41e-4  # Fluorescence lifetime [s]
        medium.get("claddingNumber").value = 21  # Physical surface tag of "CladdingShell" in the gmsh file.
        medium.get("claddingAbsorption").value = 5.5  # [1/cm]

        cross_sections_data = SpectralDecomposition(
            wavelengthsAbsorption=[900.0, 910.0],
            crossSectionAbsorption=[1.1e-21, 1.2e-21],
            wavelengthsEmission=[1020.0, 1030.0],
            crossSectionEmission=[2.0e-20, 2.48e-20],
            resolution=2,
        )
        pump = PumpProperties(
            spectralProperties=cross_sections_data,
            intensity=16e3,  # [W/cm^2]
            pumpSubsteps=100,
            wavelength=940e-9,  # [m]
            radiusX=1.5,
            radiusY=1.5,
            superGaussianOrder=40,
            myCustomVar=6,
        )

        phi_ase = PhiASE(
            spectralProperties=cross_sections_data,
            minRaysPerSample=1000,
            maxRaysPerSample=1000,
            repetitions=1,
            adaptiveSteps=1,
            mseThreshold=0.005,
            useReflections=True,
            backend="Host_Cpu_CpuSerial",
            parallelMode="single",
            numDevices=1,
        )

        simulation = Simulation(
            gainMedium=medium,
            pump=pump,
            phiASE=phi_ase,
            timeIntegrationSolver=RungeKutta4(),
            timeStep=1e-5,
            endTime=1e-3,
        )
        simulation.onInit(initFunc)
        simulation.onStep(printState)
        simulation.onStep(writeVtkState, "gmsh_minimal_phi_ase_{step:03d}.vtk")
        simulation.runSteps(3)

        last_state = simulation.getLastState()
        print(f"last completed step: {last_state.step}")


if __name__ == "__main__":
    main()
