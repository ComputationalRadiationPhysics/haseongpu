# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np

from _source_tree_import import ensure_hase_importable

ensure_hase_importable()

from HASEonGPU import (
    GainMedium,
    Grid,
    MeshTopology,
    PhiASE,
    PrimitiveFieldSpec,
    PrismSchema,
    PlanarPumpRelay,
    MonteCarloPumpSolver,
    Pump,
    PumpSpectrum,
    SuperGaussianPumpProfile,
    SurfacePumpInjector,
    RungeKutta4,
    Simulation,
    SpectralDecomposition,
    vtkWedge,
)
def initFunc(simulation):
    medium=simulation.gainMedium
def printState(state):
    print(
        f"step={state.step:03d} "
        f"time={state.time:.3e}s "
        f"mean_beta={state.beta_cells.mean():.6e} "
        f"mean_phi={state.phi_ase.mean():.6e}"
    )


def writeVtkState(state, outputFile):
    vtkWedge(outputFile, state)

def main():
    # docs:start: topology
    topology = MeshTopology.fromGrid(
        Grid(xExtent=4, yExtent=4, zExtent=0.7, tileSizeX=0.25, tileSizeZ=0.7 / 9.0)
    )
    # docs:end: topology
    # docs:start: gain-medium
    medium = GainMedium(topology=topology)
    print("betaCells shape:", medium.get("betaCells").expectedShape)

    for point in medium.getPoints():
        point.betaCells = 0.0

    for prism in medium.getPrisms():
        prism.betaVolume = 0.0

    for triangle in medium.getTriangles():
        triangle.claddingCellTypes = 0
        triangle.reflectivities = [0.0, 0.0]

    medium.get("refractiveIndices").value = np.asarray([2.0, 1.0, 3.0, 4.0], dtype=np.float32)
    medium.get("nTot").value = 1.388e20 * 2.0  # Doping density [1/cm^3]
    medium.get("crystalTFluo").value = 9.41e-4  # Fluorescence lifetime [s]
    medium.get("claddingNumber").value = 1
    medium.get("claddingAbsorption").value = 5.5  # [1/cm]

    class ThermalPrism(PrismSchema):
        temperature = PrimitiveFieldSpec(
            "temperature", "custom_temperature", np.float64, unit="K", backendRequired=False
        )

    medium.withPrimitiveSchema(ThermalPrism)
    for prism in medium.getPrisms():
        prism.temperature = 300.0

    first_prism = next(iter(medium.getPrisms()))
    print("prism fields:", first_prism.getFields())
    for field in first_prism.getFields():
        if field.name == "temperature":
            field.value(305.0)
    print("first prism temperature:", first_prism.temperature)
    # docs:end: gain-medium
    # docs:start: spectral-decomposition
    cross_sections_data = SpectralDecomposition(
        wavelengthsAbsorption=[900.0, 910.0],
        crossSectionAbsorption=[1.1e-21, 1.2e-21],
        wavelengthsEmission=[1020.0, 1030.0],
        crossSectionEmission=[2.0e-20, 2.48e-20],
        resolution=2,
    )
    print("spectral fields:", cross_sections_data.getFields())
    # docs:end: spectral-decomposition
    # docs:start: pump-properties
    pump_profile = SuperGaussianPumpProfile(radius_u=1.5, radius_v=1.5, exponent=40)
    pump = Pump(
        total_power=16e3 * 16.0,
        spectrum=PumpSpectrum.monochromatic(940e-9),
        cross_sections=cross_sections_data,
        profile=pump_profile,
    )
    pump_solver = MonteCarloPumpSolver(ray_count=100000)
    # docs:end: pump-properties


    # docs:start: phi-ase
    phi_ase = PhiASE(
        spectralProperties=cross_sections_data,
        forwardRayCount=1000,
        repetitions=1,
        relativeStandardErrorThreshold=0.1,
        useReflections=True,
        backend="Host_Cpu_CpuSerial",
        parallelMode="single",
        numDevices=1,
    )
    # docs:end: phi-ase

    # docs:start: simulation
    simulation = Simulation(
        gain_medium=medium,
        phi_ase=phi_ase,
        time_integrator=RungeKutta4(),
        time_step_size=1e-5,
        pump_solver=pump_solver,
        max_time=1e-3,
    ).add_pump(
        pump,
        injection_method=SurfacePumpInjector(surface_domains=(1,)),
        relays=(PlanarPumpRelay.retroreflect((2,)),),
    )
    simulation.on_init(initFunc)
    simulation.on_step(printState)
    simulation.on_step(writeVtkState, "minimal_phi_ase_{step:03d}.vtk")
    simulation.step(3)
    # Equivalent long run:
    # simulation.run_until(max_time=1e-3)
    # docs:end: simulation

    # docs:start: results
    last_state = simulation.get_last_state()
    print(f"last completed step: {last_state.step}")
    # docs:end: results


if __name__ == "__main__":
    main()
