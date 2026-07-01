# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

from types import SimpleNamespace

import numpy as np
import pytest

from HASEonGPU import (
    Constants,
    CrossSectionData,
    ExponentialEuler,
    ExplicitEuler,
    GainMedium,
    Heun,
    ImplicitEuler,
    Midpoint,
    PhiASE,
    PumpProperties,
    PumpRadiationProfile,
    RungeKutta4,
    Simulation,
    oneDimensionalZTraversalPumpRate,
)
from pyInclude.openpmd import transport


@pytest.fixture
def fakeCppSimulation(monkeypatch, smallTopology):
    captured = []

    def make_state(step, simulation, pump_steps):
        shape = (smallTopology.numberOfPoints, smallTopology.levels)
        volume_shape = (smallTopology.numberOfTriangles, smallTopology.levels - 1)
        pump_active = pump_steps is None or step <= pump_steps
        return SimpleNamespace(
            step=step,
            time=step * simulation.timeStep,
            betaCells=np.full(shape, 0.25 * step),
            betaVolume=np.full(volume_shape, 0.125 * step),
            phiAse=np.full(shape, float(step)),
            dndtAse=np.zeros(shape),
            dndtPump=np.ones(shape) if pump_active else np.zeros(shape),
            aseResult=object(),
        )

    def fake_run_simulation(simulation, *, steps, pumpSteps=None, transport=None):
        captured.append(
            {
                "simulation": simulation,
                "steps": steps,
                "pumpSteps": pumpSteps,
                "transport": transport,
            }
        )
        return [make_state(step, simulation, pumpSteps) for step in range(1, steps + 1)]

    monkeypatch.setattr(transport, "runSimulation", fake_run_simulation)
    return captured


def realPhiAse(crossSections, *, openpmdBackend="adios"):
    return PhiASE(spectralProperties=crossSections, openpmdBackend=openpmdBackend)


def testCompiledSimulationDelegatesRunStepsToCppTransport(
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    phi_ase = realPhiAse(crossSections)
    simulation = Simulation(
        gainMedium=smallGainMedium,
        pump=pumpProperties,
        phiASE=phi_ase,
        timeIntegrationSolver="heun",
        timeStep=1e-5,
    )

    simulation.runSteps(1, pumpSteps=0)

    state = simulation.getLastState()
    assert fakeCppSimulation == [
        {
            "simulation": simulation,
            "steps": 1,
            "pumpSteps": 0,
            "transport": "adios",
        }
    ]
    assert state.step == 1
    assert np.allclose(state.betaCells, 0.25)
    assert np.allclose(simulation.gainMedium.get("betaCells").value, 0.25)


def testTimeSteppedSimulationRunsCallbacksFromCppSnapshots(
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    seen = []
    simulation = Simulation(
        gainMedium=smallGainMedium,
        pump=pumpProperties,
        phiASE=realPhiAse(crossSections),
        timeIntegrationSolver=ExponentialEuler(),
        timeStep=1e-5,
    ).onStep(seen.append)

    simulation.runSteps(2)

    assert simulation.getLastState() is seen[-1]
    assert simulation.getResults() is seen[-1]
    assert simulation.lastState is seen[-1]
    assert len(seen) == 2
    assert seen[-1].step == 2
    assert seen[-1].time == 2e-5
    assert seen[-1].betaCells.shape == (4, 3)
    assert seen[-1].betaVolume.shape == (2, 2)


def testRunStepsCanLimitPumpContribution(
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    seen = []
    simulation = Simulation(
        gainMedium=smallGainMedium,
        pump=pumpProperties,
        phiASE=realPhiAse(crossSections),
        timeIntegrationSolver=ExponentialEuler(),
        timeStep=1e-5,
    ).onStep(seen.append)

    simulation.runSteps(3, pumpSteps=1)

    assert fakeCppSimulation[-1]["pumpSteps"] == 1
    assert np.any(seen[0].dndtPump > 0.0)
    assert np.allclose(seen[1].dndtPump, 0.0)
    assert np.allclose(seen[2].dndtPump, 0.0)


def testRunStepsUsesPumpPropertiesPumpStepsByDefault(
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    pumpProperties.customProperties["pumpSteps"] = 1
    simulation = Simulation(
        gainMedium=smallGainMedium,
        pump=pumpProperties,
        phiASE=realPhiAse(crossSections),
        timeIntegrationSolver=ExponentialEuler(),
        timeStep=1e-5,
    )

    simulation.runSteps(3)

    assert fakeCppSimulation[-1]["pumpSteps"] == 1


def testRunStepsRejectsNegativePumpSteps(
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    simulation = Simulation(
        gainMedium=smallGainMedium,
        pump=pumpProperties,
        phiASE=realPhiAse(crossSections),
        timeIntegrationSolver=ExponentialEuler(),
        timeStep=1e-5,
    )

    with pytest.raises(ValueError, match="pumpSteps"):
        simulation.runSteps(1, pumpSteps=-1)


def testOnStepPassesStateBeforeUserArguments(
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    seen = []

    def record(state, label, scale=1.0):
        seen.append((label, state.step, scale, state.phiAse.shape))

    simulation = Simulation(
        gainMedium=smallGainMedium,
        pump=pumpProperties,
        phiASE=realPhiAse(crossSections),
        timeIntegrationSolver=ExponentialEuler(),
        timeStep=1e-5,
    ).onStep(record, "vtk", scale=2.0)

    simulation.runSteps(2)

    assert seen == [("vtk", 1, 2.0, (4, 3)), ("vtk", 2, 2.0, (4, 3))]


def testInitCallbacksRunBeforeCompiledTransport(
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    events = []

    def init(simulation, label, enabled=False):
        events.append(("init", label, enabled, simulation.stepIndex))
        simulation.pump.withProperty("initialized", True)

    simulation = Simulation(
        gainMedium=smallGainMedium,
        pump=pumpProperties,
        phiASE=realPhiAse(crossSections),
        timeIntegrationSolver=ExponentialEuler(),
        timeStep=1e-5,
    ).onInit(init, "setup", enabled=True)

    simulation.runSteps(2)
    simulation.runSteps(1)

    assert events == [("init", "setup", True, 0)]
    assert simulation.pump.getProperty("initialized") is True
    assert simulation.stepIndex == 3
    assert simulation.time == 3.0000000000000004e-5


def testCompiledSimulationRejectsPythonBeforeStepCallbacks(
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    simulation = Simulation(
        gainMedium=smallGainMedium,
        pump=pumpProperties,
        phiASE=realPhiAse(crossSections),
        timeIntegrationSolver="explicit-euler",
        timeStep=1e-5,
    ).beforeStep(lambda simulation: None)

    with pytest.raises(ValueError, match="beforeStep"):
        simulation.runSteps(1)


def testCompiledSimulationRejectsExternalOpenPmdSessionOwnership(
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    simulation = Simulation(
        gainMedium=smallGainMedium,
        pump=pumpProperties,
        phiASE=realPhiAse(crossSections),
        timeIntegrationSolver=ExponentialEuler(),
        timeStep=1e-5,
    )

    with pytest.raises(ValueError, match=r"owns its C\+\+ openPMD lifetime"):
        simulation.runSteps(1, openpmdSession="persistent")


def testCompiledSimulationRejectsCustomPythonPumpSolver(
    monkeypatch,
    smallGainMedium,
    crossSections,
):
    class ConstantPumpSolver:
        def step(self, input, pump):
            return input["betaCell"] + 0.25

    monkeypatch.setattr(transport, "_ensure_backend_available", lambda backend: None)
    monkeypatch.setattr(transport, "findCalcPhiAse", lambda: "calcPhiASE")
    pump = PumpProperties(
        spectralProperties=crossSections,
        intensity=16e3,
        pumpSubsteps=100,
        solver=ConstantPumpSolver(),
        wavelength=940e-9,
    )

    simulation = Simulation(
        gainMedium=smallGainMedium,
        pump=pump,
        phiASE=realPhiAse(crossSections),
        timeIntegrationSolver=ExponentialEuler(),
        timeStep=1e-5,
    )

    with pytest.raises(ValueError, match="custom Python pump solvers"):
        simulation.step()


def testDefaultPumpSolverRequiresGaussianRadius(crossSections):
    try:
        PumpProperties(
            spectralProperties=crossSections,
            intensity=16e3,
            wavelength=940e-9,
        )
    except ValueError as exc:
        assert "radiusX" in str(exc)
    else:
        raise AssertionError("default Gaussian pump accepted missing radiusX")


def testPumpPropertiesAcceptsArbitraryDirectKeywords(crossSections):
    pump = PumpProperties(
        spectralProperties=crossSections,
        intensity=16e3,
        pumpSubsteps=100,
        wavelength=940e-9,
        radiusX=1.5,
        radiusY=1.5,
        superGaussianOrder=40,
        ji=3,
    )

    assert pump.getProperty("ji") == 3
    assert pump.ji == 3
    assert pump.radiusX == 1.5
    assert pump.superGaussianOrder == 40


def testOneDimensionalZTraversalUsesProfileCenter(
    smallTopology,
):
    spectra = CrossSectionData.monochromatic(
        wavelength=940e-9,
        crossSectionAbsorption=1.0e-22,
        crossSectionEmission=0.0,
    )
    medium = GainMedium(topology=smallTopology).withPhysicalProperties(
        betaCells=np.zeros((smallTopology.numberOfPoints, smallTopology.levels)),
        claddingCellTypes=np.zeros(smallTopology.numberOfTriangles, dtype=np.uint32),
        refractiveIndices=[1.8, 1.0, 1.8, 1.0],
        reflectivities=np.zeros((smallTopology.numberOfTriangles, 2)),
        nTot=0.0,
        crystalTFluo=9.5e-4,
        claddingNumber=1,
        claddingAbsorption=0.0,
    )
    profile = PumpRadiationProfile(
        intensity=1.0,
        wavelengths=[940e-9],
        waist=(0.2, 0.2),
        center=smallTopology.points[0],
        propagationDirection=(0.0, 0.0, 1.0),
        backReflection=False,
        superGaussianOrder=2,
    )
    pump = PumpProperties(crossSections=spectra, profile=profile)

    rate = oneDimensionalZTraversalPumpRate(
        smallTopology.points,
        medium.get("betaCells").value,
        pump,
        medium,
    )

    assert np.all(rate[0] >= rate[1:])
    assert np.any(rate[0] > rate[1:])


def testPumpRadiationProfileUsesCrystalSpectraAndDefaultUnitWeights(
    smallTopology,
):
    spectra = CrossSectionData(
        wavelengthsAbsorption=[900e-9, 940e-9],
        crossSectionAbsorption=[1.0e-22, 2.0e-22],
        wavelengthsEmission=[900e-9, 940e-9],
        crossSectionEmission=[0.0, 0.0],
    )
    medium = GainMedium(topology=smallTopology).withPhysicalProperties(
        betaCells=np.zeros((smallTopology.numberOfPoints, smallTopology.levels)),
        claddingCellTypes=np.zeros(smallTopology.numberOfTriangles, dtype=np.uint32),
        refractiveIndices=[1.8, 1.0, 1.8, 1.0],
        reflectivities=np.zeros((smallTopology.numberOfTriangles, 2)),
        nTot=0.0,
        crystalTFluo=9.5e-4,
        claddingNumber=1,
        claddingAbsorption=0.0,
    )
    profile = PumpRadiationProfile(
        intensity=10.0,
        wavelengths=[900e-9, 940e-9],
        waist=(1.5, 1.5),
        propagationDirection=(0.0, 0.0, 1.0),
        backReflection=False,
        superGaussianOrder=40,
    )
    pump = PumpProperties(
        crossSections=spectra,
        profile=profile,
    )

    rate = oneDimensionalZTraversalPumpRate(
        smallTopology.points,
        medium.get("betaCells").value,
        pump,
        medium,
        Constants(c=3.0e8, h=6.0e-34),
    )

    points = np.asarray(smallTopology.points, dtype=np.float64)
    intensity = 10.0 * np.exp(-((points[:, 0] ** 2 + points[:, 1] ** 2) / (1.5**2)) ** 20.0)
    expected_per_point = intensity * (
        1.0e-22 * 900e-9 / (6.0e-34 * 3.0e8)
        + 2.0e-22 * 940e-9 / (6.0e-34 * 3.0e8)
    )

    assert np.allclose(rate, expected_per_point[:, None])


def testOneDimensionalZTraversalSupportsMultichromaticWeightsAndDirection(
    smallTopology,
):
    spectra = CrossSectionData(
        wavelengthsAbsorption=[900e-9, 940e-9],
        crossSectionAbsorption=[1.0e-22, 2.0e-22],
        wavelengthsEmission=[900e-9, 940e-9],
        crossSectionEmission=[0.5e-22, 1.0e-22],
    )
    medium = GainMedium(topology=smallTopology).withPhysicalProperties(
        betaCells=np.full((smallTopology.numberOfPoints, smallTopology.levels), 0.25),
        claddingCellTypes=np.zeros(smallTopology.numberOfTriangles, dtype=np.uint32),
        refractiveIndices=[1.8, 1.0, 1.8, 1.0],
        reflectivities=np.zeros((smallTopology.numberOfTriangles, 2)),
        nTot=0.0,
        crystalTFluo=9.5e-4,
        claddingNumber=1,
        claddingAbsorption=0.0,
    )
    pump = PumpProperties(
        spectralProperties=spectra,
        intensity=10.0,
        radiusX=1.5,
        backReflection=True,
        reflectivity=0.5,
        propagationDirection=(0.0, 0.0, -1.0),
        spectralWeights=[0.25, 0.75],
    )

    rate = oneDimensionalZTraversalPumpRate(
        smallTopology.points,
        medium.get("betaCells").value,
        pump,
        medium,
        Constants(c=3.0e8, h=6.0e-34),
    )

    points = np.asarray(smallTopology.points, dtype=np.float64)
    intensity = 10.0 * np.exp(-((points[:, 0] ** 2 + points[:, 1] ** 2) / (1.5**2)) ** 20.0)
    wavelengths = np.asarray([900e-9, 940e-9])
    sigma_abs = np.asarray([1.0e-22, 2.0e-22])
    sigma_ems = np.asarray([0.5e-22, 1.0e-22])
    weights = np.asarray([0.25, 0.75])
    beta = 0.25
    expected_per_point = np.sum(
        weights
        * (sigma_abs - beta * (sigma_abs + sigma_ems))
        * (1.0 + 0.5)
        * intensity[:, None]
        * wavelengths
        / (6.0e-34 * 3.0e8),
        axis=1,
    )

    assert np.allclose(rate, expected_per_point[:, None])


def testTimeIntegrationSolverIsMandatory(
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    try:
        Simulation(
            gainMedium=smallGainMedium,
            pump=pumpProperties,
            phiASE=realPhiAse(crossSections),
            timeIntegrationSolver=None,
            timeStep=1e-5,
        )
    except ValueError as exc:
        assert "timeIntegrationSolver" in str(exc)
    else:
        raise AssertionError("Simulation accepted a missing timeIntegrationSolver")


def testTimeIntegrationSolversCanStepSimulation(
    fakeCppSimulation,
    pumpProperties,
    crossSections,
    smallTopology,
):
    solvers = [ExplicitEuler(), Heun(), Midpoint(), RungeKutta4(), ImplicitEuler(iterations=2)]

    for solver in solvers:
        medium = GainMedium(topology=smallTopology).withPhysicalProperties(
            betaCells=np.zeros((smallTopology.numberOfPoints, smallTopology.levels)),
            claddingCellTypes=np.zeros(smallTopology.numberOfTriangles, dtype=np.uint32),
            refractiveIndices=[1.8, 1.0, 1.8, 1.0],
            reflectivities=np.zeros((smallTopology.numberOfTriangles, 2)),
            nTot=2.76e20,
            crystalTFluo=9.5e-4,
            claddingNumber=1,
            claddingAbsorption=0.0,
        )
        state = Simulation(
            gainMedium=medium,
            pump=pumpProperties,
            phiASE=realPhiAse(crossSections),
            timeIntegrationSolver=solver,
            timeStep=1e-5,
        ).step()

        assert state.betaCells.shape == (smallTopology.numberOfPoints, smallTopology.levels)
        assert np.all(np.isfinite(state.betaCells))
