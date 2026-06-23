# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later


import numpy as np

from HASEonGPU import (
    Constants,
    ExponentialEuler,
    ExplicitEuler,
    GainMedium,
    Heun,
    ImplicitEuler,
    Midpoint,
    PumpProperties,
    RungeKutta4,
    Simulation,
)


def testTimeSteppedSimulationRunsCallbacksWithFakeAse(
    smallGainMedium,
    smallTopology,
    pumpProperties,
    makeFakePhiAse,
):
    phiAse = makeFakePhiAse(smallTopology)
    seen = []
    simulation = Simulation(
        gainMedium=smallGainMedium,
        pump=pumpProperties,
        phiASE=phiAse,
        timeIntegrationSolver=ExponentialEuler(),
        timeStep=1e-5,
    ).onStep(seen.append)

    simulation.runSteps(2)

    assert simulation.getLastState() is seen[-1]
    assert simulation.getResults() is seen[-1]
    assert simulation.lastState is seen[-1]
    assert len(seen) == 2
    assert seen[-1].step == 2
    assert seen[-1].betaCells.shape == (4, 3)
    assert seen[-1].betaVolume.shape == (2, 2)


def testRunStepsCanLimitPumpContribution(
    smallGainMedium,
    smallTopology,
    pumpProperties,
    makeFakePhiAse,
):
    class ConstantPumpSolver:
        def step(self, input, pump):
            return np.asarray(input["betaCell"], dtype=np.float64) + 1.0e-6

    pumpProperties.customProperties["solver"] = ConstantPumpSolver()
    seen = []
    simulation = Simulation(
        gainMedium=smallGainMedium,
        pump=pumpProperties,
        phiASE=makeFakePhiAse(smallTopology),
        timeIntegrationSolver=ExponentialEuler(),
        timeStep=1e-5,
    ).onStep(seen.append)

    simulation.runSteps(3, pumpSteps=1)

    assert np.any(seen[0].dndtPump > 0.0)
    assert np.allclose(seen[1].dndtPump, 0.0)
    assert np.allclose(seen[2].dndtPump, 0.0)


def testRunStepsUsesPumpPropertiesPumpStepsByDefault(
    smallGainMedium,
    smallTopology,
    pumpProperties,
    makeFakePhiAse,
):
    class ConstantPumpSolver:
        def step(self, input, pump):
            return np.asarray(input["betaCell"], dtype=np.float64) + 1.0e-6

    pumpProperties.customProperties["solver"] = ConstantPumpSolver()
    pumpProperties.customProperties["pumpSteps"] = 1
    seen = []
    simulation = Simulation(
        gainMedium=smallGainMedium,
        pump=pumpProperties,
        phiASE=makeFakePhiAse(smallTopology),
        timeIntegrationSolver=ExponentialEuler(),
        timeStep=1e-5,
    ).onStep(seen.append)

    simulation.runSteps(3)

    assert np.any(seen[0].dndtPump > 0.0)
    assert np.allclose(seen[1].dndtPump, 0.0)
    assert np.allclose(seen[2].dndtPump, 0.0)


def testOnStepPassesStateBeforeUserArguments(
    smallGainMedium,
    smallTopology,
    pumpProperties,
    makeFakePhiAse,
):
    seen = []

    def record(state, label, scale=1.0):
        seen.append((label, state.step, scale, state.phiAse.shape))

    simulation = Simulation(
        gainMedium=smallGainMedium,
        pump=pumpProperties,
        phiASE=makeFakePhiAse(smallTopology),
        timeIntegrationSolver=ExponentialEuler(),
        timeStep=1e-5,
    ).onStep(record, "vtk", scale=2.0)

    simulation.runSteps(2)

    assert seen == [("vtk", 1, 2.0, (4, 3)), ("vtk", 2, 2.0, (4, 3))]


def testInitAndBeforeStepPassSimulationBeforeUserArguments(
    smallGainMedium,
    smallTopology,
    pumpProperties,
    makeFakePhiAse,
):
    events = []

    def init(simulation, label, enabled=False):
        events.append(("init", label, enabled, simulation.stepIndex))

    def before(simulation, label, scale=1.0):
        events.append(("before", label, scale, simulation.stepIndex))

    simulation = Simulation(
        gainMedium=smallGainMedium,
        pump=pumpProperties,
        phiASE=makeFakePhiAse(smallTopology),
        timeIntegrationSolver=ExponentialEuler(),
        timeStep=1e-5,
    ).onInit(init, "setup", enabled=True).beforeStep(before, "pre", scale=3.0)

    simulation.runSteps(2)

    assert events == [
        ("init", "setup", True, 0),
        ("before", "pre", 3.0, 0),
        ("before", "pre", 3.0, 1),
    ]


def testTimeSteppedSimulationLifecycleHooksGetSimulation(
    smallGainMedium,
    smallTopology,
    pumpProperties,
    makeFakePhiAse,
):
    phiAse = makeFakePhiAse(smallTopology)
    events = []

    def init(simulation):
        events.append(("init", simulation.stepIndex, simulation.time, simulation.gainMedium is smallGainMedium))
        simulation.pump.withProperty("initialized", True)

    def beforeStep(simulation):
        events.append(("before", simulation.stepIndex, simulation.time, simulation.pump.getProperty("initialized")))

    def afterStep(state):
        events.append(("after", state.step, state.time))

    simulation = Simulation(
        gainMedium=smallGainMedium,
        pump=pumpProperties,
        phiASE=phiAse,
        timeIntegrationSolver=ExponentialEuler(),
        timeStep=1e-5,
    ).onInit(init).beforeStep(beforeStep).onStep(afterStep)

    simulation.runSteps(2)
    simulation.runSteps(1)

    assert events == [
        ("init", 0, 0.0, True),
        ("before", 0, 0.0, True),
        ("after", 1, 1e-5),
        ("before", 1, 1e-5, True),
        ("after", 2, 2e-5),
        ("before", 2, 2e-5, True),
        ("after", 3, 3.0000000000000004e-5),
    ]


def testTimeSteppedSimulationAcceptsCustomPumpSolver(
    smallGainMedium,
    smallTopology,
    crossSections,
    makeFakePhiAse,
):
    class ConstantPumpSolver:
        def step(self, input, pump):
            return input["betaCell"] + 0.25

    phiAse = makeFakePhiAse(smallTopology)
    pump = PumpProperties(
        spectralProperties=crossSections,
        intensity=16e3,
        pumpSubsteps=100,
        solver=ConstantPumpSolver(),
        wavelength=940e-9,
    )

    state = Simulation(
        gainMedium=smallGainMedium,
        pump=pump,
        phiASE=phiAse,
        timeIntegrationSolver=ExponentialEuler(),
        timeStep=1e-5,
    ).step()

    tau = 9.5e-4
    expected = tau * (0.25 / 1e-5) * (1.0 - np.exp(-1e-5 / tau))
    assert np.allclose(state.betaCells, expected)


def testCustomPumpSolverCanReadAdditionalPumpProperties(
    smallGainMedium,
    smallTopology,
    crossSections,
    makeFakePhiAse,
):
    class ScaledPumpSolver:
        def step(self, input, pump):
            return input["betaCell"] + pump.getProperty("betaIncrement")

    phiAse = makeFakePhiAse(smallTopology)
    pump = PumpProperties(
        spectralProperties=crossSections,
        intensity=16e3,
        pumpSubsteps=100,
        solver=ScaledPumpSolver(),
        wavelength=940e-9,
        customProperties={"betaIncrement": 0.1},
    ).withProperty("profileName", "user-defined")

    state = Simulation(
        gainMedium=smallGainMedium,
        pump=pump,
        phiASE=phiAse,
        timeIntegrationSolver=ExponentialEuler(),
        timeStep=1e-5,
    ).step()

    tau = 9.5e-4
    expected = tau * (0.1 / 1e-5) * (1.0 - np.exp(-1e-5 / tau))
    assert pump.getProperty("profileName") == "user-defined"
    assert np.allclose(state.betaCells, expected)


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


def testAseDerivativeUsesPhiAseScalingFromBackend(
    smallGainMedium,
    smallTopology,
    pumpProperties,
    crossSections,
    makeFakePhiAse,
):
    simulation = Simulation(
        gainMedium=smallGainMedium,
        pump=pumpProperties,
        phiASE=makeFakePhiAse(smallTopology),
        timeIntegrationSolver=ExponentialEuler(),
        timeStep=1e-5,
        crossSections=crossSections,
    )
    betaCells = np.full((smallTopology.numberOfPoints, smallTopology.levels), 0.25)
    phiAse = np.full_like(betaCells, 2.0)

    derivative = simulation._aseDerivative(phiAse, betaCells=betaCells)
    gainPerDensity = 0.25 * (
        crossSections.emissionAt(940e-9) + crossSections.absorptionAt(940e-9)
    ) - crossSections.absorptionAt(940e-9)

    assert np.allclose(derivative, gainPerDensity * phiAse)


def testPumpPropertiesAcceptsArbitraryDirectKeywords(crossSections):
    pump = PumpProperties(
        spectralProperties=crossSections,
        intensity=16e3,
        pumpSubsteps=100,
        wavelength=940e-9,
        radiusX=1.5,
        radiusY=1.5,
        exponent=40,
        ji=3,
    )

    assert pump.getProperty("ji") == 3
    assert pump.ji == 3
    assert pump.radiusX == 1.5
    assert pump.exponent == 40


def testTimeIntegrationSolverIsMandatory(
    smallGainMedium,
    pumpProperties,
    makeFakePhiAse,
    smallTopology,
):
    phiAse = makeFakePhiAse(smallTopology)

    try:
        Simulation(
            gainMedium=smallGainMedium,
            pump=pumpProperties,
            phiASE=phiAse,
            timeIntegrationSolver=None,
            timeStep=1e-5,
        )
    except ValueError as exc:
        assert "timeIntegrationSolver" in str(exc)
    else:
        raise AssertionError("Simulation accepted a missing timeIntegrationSolver")


def testTimeIntegrationSolversCanStepSimulation(
    smallGainMedium,
    pumpProperties,
    makeFakePhiAse,
    smallTopology,
):
    solvers = [ExplicitEuler(), Heun(), Midpoint(), RungeKutta4(), ImplicitEuler(iterations=2)]

    for solver in solvers:
        medium = GainMedium(topology=smallTopology).withPhysicalProperties(
            betaCells=np.zeros((smallTopology.numberOfPoints, smallTopology.levels)),
            claddingCellTypes=np.zeros(smallTopology.numberOfTriangles, dtype=np.uint32),
            refractiveIndices=[1.8, 1.0, 1.8, 1.0],
            reflectivities=np.zeros((2, smallTopology.numberOfTriangles)),
            nTot=2.76e20,
            crystalTFluo=9.5e-4,
            claddingNumber=1,
            claddingAbsorption=0.0,
        )
        state = Simulation(
            gainMedium=medium,
            pump=pumpProperties,
            phiASE=makeFakePhiAse(smallTopology),
            timeIntegrationSolver=solver,
            timeStep=1e-5,
        ).step()

        assert state.betaCells.shape == (smallTopology.numberOfPoints, smallTopology.levels)
        assert np.all(np.isfinite(state.betaCells))
