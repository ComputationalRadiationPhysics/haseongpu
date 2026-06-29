# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later


import numpy as np
import pytest

from HASEonGPU import (
    Constants,
    OneDimensionalZTraversal,
    CrossSectionData,
    ExponentialEuler,
    ExplicitEuler,
    FrozenPhiAseRungeKutta4,
    GainMedium,
    Heun,
    ImplicitEuler,
    Midpoint,
    PumpProperties,
    RungeKutta4,
    Simulation,
    PumpRadiationProfile,
    oneDimensionalZTraversalPumpRate,
)


@pytest.fixture(autouse=True)
def _use_file_openpmd_backend(monkeypatch):
    monkeypatch.setenv("HASE_OPENPMD_BACKEND", "adios")


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


def testSimulationCanDisableAseWithoutBreakingState(
    smallGainMedium,
    smallTopology,
    pumpProperties,
    makeFakePhiAse,
):
    phiAse = makeFakePhiAse(smallTopology)

    def failRun(*args, **kwargs):
        raise AssertionError("PhiASE.run should not be called when enableAse=False")

    phiAse.run = failRun
    state = Simulation(
        gainMedium=smallGainMedium,
        pump=pumpProperties,
        phiASE=phiAse,
        timeIntegrationSolver=ExponentialEuler(),
        timeStep=1e-5,
        enableAse=False,
    ).step()

    assert state.aseResult is None
    assert state.phiAse.shape == (smallTopology.numberOfPoints, smallTopology.levels)
    assert state.dndtAse.shape == (smallTopology.numberOfPoints, smallTopology.levels)
    assert np.allclose(state.phiAse, 0.0)
    assert np.allclose(state.dndtAse, 0.0)
    assert np.all(np.isfinite(state.betaCells))


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
    assert np.allclose(state.betaCells[:, -1], expected)


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


def testOneDimensionalZTraversalProducesFrozenStatePumpRate(
    smallGainMedium,
    smallTopology,
    crossSections,
    makeFakePhiAse,
):
    pump = PumpProperties(
        spectralProperties=crossSections,
        intensity=16e3,
        wavelength=940e-9,
        radiusX=1.5,
        backReflection=False,
        propagationDirection=(0.0, 0.0, 1.0),
        solver=OneDimensionalZTraversal(),
    )
    state = Simulation(
        gainMedium=smallGainMedium,
        pump=pump,
        phiASE=makeFakePhiAse(smallTopology),
        timeIntegrationSolver=ExplicitEuler(),
        timeStep=1e-5,
        enableAse=False,
    ).step()

    expected = oneDimensionalZTraversalPumpRate(
        smallTopology.points,
        np.zeros((smallTopology.numberOfPoints, smallTopology.levels)),
        pump,
        smallGainMedium,
    )

    assert np.allclose(state.dndtPump, expected)
    assert np.all(state.dndtPump >= 0.0)
    assert np.allclose(state.phiAse, 0.0)


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
        reflectivities=np.zeros((2, smallTopology.numberOfTriangles)),
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
    pump = PumpProperties(crossSections=spectra, profile=profile, solver=OneDimensionalZTraversal())

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
        reflectivities=np.zeros((2, smallTopology.numberOfTriangles)),
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
        solver=OneDimensionalZTraversal(),
    )

    rate = oneDimensionalZTraversalPumpRate(
        smallTopology.points,
        medium.get("betaCells").value,
        pump,
        medium,
        Constants(c=3.0e8, h=6.0e-34),
    )

    points = np.asarray(smallTopology.points, dtype=np.float64)
    intensity = 10.0 * np.exp(-((points[:, 0] ** 2 + points[:, 1] ** 2) / (1.5 ** 2)) ** 20.0)
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
        reflectivities=np.zeros((2, smallTopology.numberOfTriangles)),
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
        solver=OneDimensionalZTraversal(),
    )

    rate = oneDimensionalZTraversalPumpRate(
        smallTopology.points,
        medium.get("betaCells").value,
        pump,
        medium,
        Constants(c=3.0e8, h=6.0e-34),
    )

    points = np.asarray(smallTopology.points, dtype=np.float64)
    intensity = 10.0 * np.exp(-((points[:, 0] ** 2 + points[:, 1] ** 2) / (1.5 ** 2)) ** 20.0)
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
        superGaussianOrder=40,
        ji=3,
    )

    assert pump.getProperty("ji") == 3
    assert pump.ji == 3
    assert pump.radiusX == 1.5
    assert pump.superGaussianOrder == 40


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


def _mediumLikeSmallTopology(smallTopology, betaCells=None):
    if betaCells is None:
        betaCells = np.zeros((smallTopology.numberOfPoints, smallTopology.levels))
    return GainMedium(topology=smallTopology).withPhysicalProperties(
        betaCells=betaCells,
        claddingCellTypes=np.zeros(smallTopology.numberOfTriangles, dtype=np.uint32),
        refractiveIndices=[1.8, 1.0, 1.8, 1.0],
        reflectivities=np.zeros((2, smallTopology.numberOfTriangles)),
        nTot=2.76e20,
        crystalTFluo=9.5e-4,
        claddingNumber=1,
        claddingAbsorption=0.0,
    )


def testFrozenPhiAseRungeKutta4RunsAseBackendOncePerStep(
    smallTopology,
    pumpProperties,
    makeFakePhiAse,
):
    phiAse = makeFakePhiAse(smallTopology)
    calls = []

    def run(*args, **kwargs):
        calls.append(np.asarray(kwargs["gainMedium"].get("betaCells").value).copy())
        return phiAse

    phiAse.run = run
    state = Simulation(
        gainMedium=_mediumLikeSmallTopology(smallTopology),
        pump=pumpProperties,
        phiASE=phiAse,
        timeIntegrationSolver=FrozenPhiAseRungeKutta4(),
        timeStep=1e-5,
    ).step()

    assert len(calls) == 1
    assert state.aseResult is not None
    assert state.phiAse.shape == (smallTopology.numberOfPoints, smallTopology.levels)


def testSimulationPrePumpSkipsAseBackendForInitialStepOnly(
    smallTopology,
    pumpProperties,
    makeFakePhiAse,
):
    class ConstantPumpSolver:
        def step(self, input, pump):
            return np.asarray(input["betaCell"], dtype=np.float64) + 1.0e-6

    pumpProperties.customProperties["solver"] = ConstantPumpSolver()
    phiAse = makeFakePhiAse(smallTopology)
    calls = []

    def run(*args, **kwargs):
        calls.append(np.asarray(kwargs["gainMedium"].get("betaCells").value).copy())
        return phiAse

    phiAse.run = run
    simulation = Simulation(
        gainMedium=_mediumLikeSmallTopology(smallTopology),
        pump=pumpProperties,
        phiASE=phiAse,
        timeIntegrationSolver=FrozenPhiAseRungeKutta4(),
        timeStep=1e-5,
        enableAse=True,
        prePump=True,
    )

    first = simulation.step()
    second = simulation.step()

    assert simulation.enableAse is True
    assert first.aseResult is None
    assert np.allclose(first.phiAse, 0.0)
    assert np.allclose(first.dndtAse, 0.0)
    assert second.aseResult is not None
    assert len(calls) == 1
    assert np.any(calls[0] > 0.0)


def testFrozenPhiAseRungeKutta4SkipsAseBackendWhenDisabled(
    smallTopology,
    pumpProperties,
    makeFakePhiAse,
):
    class LinearPumpSolver:
        def step(self, input, pump):
            beta = np.asarray(input["betaCell"], dtype=np.float64)
            return beta + input["_timeStep"] * (0.2 + 0.05 * beta)

    pumpProperties.customProperties["solver"] = LinearPumpSolver()
    initial = np.full((smallTopology.numberOfPoints, smallTopology.levels), 0.1)

    frozenPhiAse = makeFakePhiAse(smallTopology)
    frozenPhiAse.run = lambda *args, **kwargs: (_ for _ in ()).throw(
        AssertionError("PhiASE.run should not be called when enableAse=False")
    )
    frozenState = Simulation(
        gainMedium=_mediumLikeSmallTopology(smallTopology, betaCells=initial.copy()),
        pump=pumpProperties,
        phiASE=frozenPhiAse,
        timeIntegrationSolver=FrozenPhiAseRungeKutta4(),
        timeStep=1e-5,
        enableAse=False,
    ).step()

    regularPhiAse = makeFakePhiAse(smallTopology)
    regularPhiAse.run = lambda *args, **kwargs: (_ for _ in ()).throw(
        AssertionError("PhiASE.run should not be called when enableAse=False")
    )
    regularState = Simulation(
        gainMedium=_mediumLikeSmallTopology(smallTopology, betaCells=initial.copy()),
        pump=pumpProperties,
        phiASE=regularPhiAse,
        timeIntegrationSolver=RungeKutta4(),
        timeStep=1e-5,
        enableAse=False,
    ).step()

    assert np.allclose(frozenState.betaCells, regularState.betaCells)
    assert np.allclose(frozenState.phiAse, 0.0)
    assert np.allclose(frozenState.dndtAse, 0.0)


def testFrozenPhiAseDerivativeReusesPhiAseWhileAseCouplingFollowsBeta(
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
        timeIntegrationSolver=FrozenPhiAseRungeKutta4(),
        timeStep=1e-5,
        crossSections=crossSections,
    )
    phiAse = np.full((smallTopology.numberOfPoints, smallTopology.levels), 2.0)
    lowBeta = np.full_like(phiAse, 0.1)
    highBeta = np.full_like(phiAse, 0.4)

    low = simulation._timeDerivativeWithFrozenPhiAse(lowBeta, 0.0, phiAse, object())
    high = simulation._timeDerivativeWithFrozenPhiAse(highBeta, 0.0, phiAse, object())

    assert np.allclose(low.phiAse, phiAse)
    assert np.allclose(high.phiAse, phiAse)
    assert not np.array_equal(low.dndtAse, high.dndtAse)


def testTimeIntegrationSolversCanStepSimulation(
    smallGainMedium,
    pumpProperties,
    makeFakePhiAse,
    smallTopology,
):
    solvers = [
        ExplicitEuler(),
        Heun(),
        Midpoint(),
        RungeKutta4(),
        FrozenPhiAseRungeKutta4(),
        ImplicitEuler(iterations=2),
    ]

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
