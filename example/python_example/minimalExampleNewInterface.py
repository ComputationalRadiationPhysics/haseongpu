# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np

from HASEonGPU import (
    GainMedium,
    Grid,
    MeshTopology,
    PhiASE,
    PumpProperties,
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
        f"mean_beta={state.betaCells.mean():.6e} "
        f"mean_phi={state.phiAse.mean():.6e}"
    )


def writeVtkState(state, outputFile):
    vtkWedge(outputFile, state)
# docs:start: custom-pump-solver
class MyPumpSolver:
    def step(self, input, pump):
        beta = input["betaCell"]
        mycustom = pump.getProperty("myCustomVar")
        pump.withProperty("myCustomVar",mycustom+1)
        return np.ones_like(beta) - beta
# docs:end: custom-pump-solver

def main():
    # docs:start: topology
    topology = MeshTopology.fromGrid(
        Grid(xExtent=4, yExtent=4, zExtent=0.7, tileSizeX=0.25, tileSizeZ=0.7 / 9.0)
    )
    # docs:end: topology
    # docs:start: gain-medium
    medium = GainMedium(topology=topology)
    print(medium.get("betaCells").expectedShape)
    medium.withPhysicalProperties(
        betaCells=np.zeros(medium.get("betaCells").expectedShape),
        claddingCellTypes=np.zeros(medium.get("claddingCellTypes").expectedShape, dtype=np.uint32),
        refractiveIndices=[2, 1, 3, 4],
        reflectivities=np.zeros(medium.get("reflectivities").expectedShape),
        nTot=1.388e20 * 2.0, # Dotierungsdichte [1/cm^3]
        crystalTFluo=9.41e-4, # fluorescence lifetime
        claddingNumber=1,
        claddingAbsorption=5.5, # [1/cm]
    )
    # docs:end: gain-medium
    # docs:start: spectral-decomposition
    cross_sections_data = SpectralDecomposition(
        wavelengthsAbsorption=[900.0, 910.0],
        crossSectionAbsorption=[1.1e-21, 1.2e-21],
        wavelengthsEmission=[1020.0, 1030.0],
        crossSectionEmission=[2.0e-20, 2.48e-20],
        resolution=2,
    )
    # docs:end: spectral-decomposition
    # docs:start: pump-properties
    pump = PumpProperties(
        spectralProperties=cross_sections_data,
        intensity=16e3, # [W/cm^2]
        pumpSubsteps=100,
        wavelength=940e-9, # [m]
        solver=MyPumpSolver(),
        radiusX=1.5,
        radiusY=1.5,
        exponent=40,
        myCustomVar=6
    )
    # docs:end: pump-properties


    # docs:start: phi-ase
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
    # docs:end: phi-ase

    # docs:start: simulation
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
    simulation.onStep(writeVtkState, "minimal_phi_ase_{step:03d}.vtk")
    simulation.runSteps(3)
    # Equivalent long run:
    # simulation.runUntil(endtime=1e-3)
    # docs:end: simulation

    # docs:start: results
    last_state = simulation.getLastState()
    print(f"last completed step: {last_state.step}")
    # docs:end: results


if __name__ == "__main__":
    main()
