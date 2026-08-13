# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import os
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

from HASEonGPU import (
    CrossSectionData,
    ExponentialEuler,
    ExplicitEuler,
    FrozenPhiAseRungeKutta4,
    GainMedium,
    Heun,
    ImplicitEuler,
    Midpoint,
    PhiASE,
    Pump,
    PumpSpectrum,
    RungeKutta4,
    SurfacePumpInjector,
    Simulation,
    VolumeTopology,
    autonomous_final,
)
import pyInclude.simulation as simulation_module
from pyInclude.openpmd import transport


@pytest.fixture
def fakeCppSimulation(monkeypatch, smallTopology):
    captured = []

    def make_state(step, simulation):
        volume_shape = (smallTopology.numberOfTriangles, smallTopology.levels - 1)
        simulation_step = simulation.current_step + step - 1
        pump_active = any(
            simulation_step < source.pumpSteps for source in simulation.pump.sources
        )
        fields = set(simulation.outputFields)
        return SimpleNamespace(
            step=step,
            time=step * simulation.timeStep,
            betaVolume=np.full(volume_shape, 0.125 * step) if "beta_volume" in fields else None,
            phiAse=np.full(volume_shape, float(step)) if "phi_ase" in fields else None,
            standardError=np.zeros(volume_shape) if "standard_error" in fields else None,
            relativeStandardError=np.zeros(volume_shape) if "relative_standard_error" in fields else None,
            totalRays=np.zeros(volume_shape, dtype=np.uint32) if "total_rays" in fields else None,
            dndtAse=np.zeros(volume_shape) if "dndt_ase" in fields else None,
            dndtPump=(np.ones(volume_shape) if pump_active else np.zeros(volume_shape))
            if "dndt_pump" in fields else None,
            aseResult=object(),
        )

    def fake_run_simulation(
        simulation,
        *,
        steps,
        transport=None,
        command_prefix=None,
        workspace_dir=None,
        on_state=None,
    ):
        call = {
            "simulation": simulation,
            "steps": steps,
            "transport": transport,
        }
        if command_prefix is not None:
            call["command_prefix"] = command_prefix
        if workspace_dir is not None:
            call["workspace_dir"] = workspace_dir
        captured.append(call)
        emitted_steps = (
            range(1, steps + 1)
            if simulation.executionMode == "synchronized-debug" or simulation.outputSteps is None
            else simulation.outputSteps
        )
        states = [make_state(step, simulation) for step in emitted_steps]
        if on_state is not None:
            for state in states:
                on_state(state)
        return states

    monkeypatch.setattr(transport, "runSimulation", fake_run_simulation)
    return captured


def realPhiAse(crossSections, *, openpmdBackend="adios"):
    return PhiASE(spectralProperties=crossSections, openpmdBackend=openpmdBackend, ase_steps=100)


def configuredSimulation(pumpSetup, **kwargs):
    return Simulation(**kwargs).add_pump(
        pumpSetup.physical, injection_method=pumpSetup.injector
    )


def testCompiledSimulationDelegatesRunStepsToCppTransport(
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    phi_ase = realPhiAse(crossSections)
    simulation = configuredSimulation(pumpProperties,
        gain_medium=smallGainMedium,
        phi_ase=phi_ase,
        time_integrator="heun",
        time_step_size=1e-5,
    )

    simulation.runSteps(1)

    state = simulation.getLastState()
    assert fakeCppSimulation == [
        {
            "simulation": simulation,
            "steps": 1,
            "transport": "adios",
        }
    ]
    assert state.step == 1
    assert np.allclose(state.betaVolume, 0.125)
    assert np.allclose(simulation.gainMedium.get("betaVolume").value, 0.125)


def testCompiledSimulationRejectsUnavailableAlpakaBackendBeforeTransport(
    monkeypatch,
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    monkeypatch.setattr(simulation_module.AlpakaBackends, "all", lambda: ["Host_Cpu_CpuSerial"])
    phi_ase = realPhiAse(crossSections)
    phi_ase.backend = "Host_Cpu_CpuOmpBlocks"
    simulation = configuredSimulation(
        pumpProperties,
        gain_medium=smallGainMedium,
        phi_ase=phi_ase,
        time_integrator="heun",
        time_step_size=1e-5,
    )

    with pytest.raises(RuntimeError, match="Host_Cpu_CpuOmpBlocks"):
        simulation.step()

    assert fakeCppSimulation == []


def testCompiledSimulationRejectsUnavailableOpenPmdBackendBeforeTransport(
    monkeypatch,
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    monkeypatch.setattr(simulation_module.AlpakaBackends, "all", lambda: ["Host_Cpu_CpuSerial"])

    def reject_openpmd(backend):
        raise RuntimeError(f"openPMD backend '{backend}' is unavailable; available backends: adios")

    monkeypatch.setattr(transport, "_ensure_backend_available", reject_openpmd)
    phi_ase = realPhiAse(crossSections, openpmdBackend="hdf5")
    phi_ase.backend = "Host_Cpu_CpuSerial"
    simulation = configuredSimulation(
        pumpProperties,
        gain_medium=smallGainMedium,
        phi_ase=phi_ase,
        time_integrator="heun",
        time_step_size=1e-5,
    )

    with pytest.raises(RuntimeError, match="available backends: adios"):
        simulation.step()

    assert fakeCppSimulation == []


def testCompiledSimulationUsesPhiAseMpiLaunchOptions(
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
    monkeypatch,
    tmp_path,
):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("HASE_MPIEXEC_EXTRA_ARGS", "--oversubscribe")
    phi_ase = PhiASE(
        spectralProperties=crossSections,
        openpmdBackend="adios",
        parallelMode="mpi",
        nPerNode=3,
    )
    simulation = configuredSimulation(
        pumpProperties,
        gain_medium=smallGainMedium,
        phi_ase=phi_ase,
        time_integrator="heun",
        time_step_size=1e-5,
    )

    simulation.runSteps(1)

    assert fakeCppSimulation[-1]["command_prefix"] == [
        "mpiexec",
        "--oversubscribe",
        "-npernode",
        "3",
    ]
    assert fakeCppSimulation[-1]["workspace_dir"] == tmp_path / "IO" / "phiase_mpi"


@pytest.mark.integration
def testCompiledSimulationMpiRanksShareOneDeviceAndAdvanceAse(
    pumpProperties,
    crossSections,
    openPmdRuntimeBackend,
    alpakaRuntimeBackend,
    monkeypatch,
    tmp_path,
):
    rank_count_text = os.environ.get("HASE_MPI_TEST_RANKS", "").strip()
    if not rank_count_text:
        pytest.skip("HASE_MPI_TEST_RANKS is not configured for this test run")
    rank_count = int(rank_count_text)

    topology = VolumeTopology.fromTetrahedra(
        np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]
        ),
        np.array([[0, 1, 2, 3]], dtype=np.uint32),
    ).withDomains(surfaceDomains={"where": "all_exterior", "domain": 1})
    gain_medium = GainMedium(topology).withPhysicalProperties(
        betaVolume=np.array([0.1]),
        claddingCellTypes=np.array([0], dtype=np.uint32),
        nTot=2.76e20,
        crystalTFluo=9.5e-4,
        claddingNumber=1,
        claddingAbsorption=0.0,
    )

    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("HASE_MPIEXEC_EXTRA_ARGS", "--oversubscribe")
    phi_ase = PhiASE(
        spectralProperties=crossSections,
        backend=alpakaRuntimeBackend,
        openpmdBackend=openPmdRuntimeBackend,
        parallelMode="mpi",
        nPerNode=rank_count,
        numDevices=1,
        minRays=256,
        maxRays=256,
        adaptiveSteps=1,
        rngSeed=1234,
        useReflections=False,
        ase_steps=2,
    )
    simulation = configuredSimulation(
        pumpProperties,
        gain_medium=gain_medium,
        phi_ase=phi_ase,
        time_integrator="explicit-euler",
        time_step_size=1e-5,
        pre_pump=True,
        output_steps=(2,),
    )

    simulation.runSteps(2)

    state = simulation.getLastState()
    assert state.step == 2
    assert np.any(np.asarray(state.phiAse) > 0.0)
    assert np.all(np.isfinite(state.betaVolume))


def testTimeSteppedSimulationRunsCallbacksFromEveryDefaultSnapshot(
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    seen = []
    simulation = configuredSimulation(pumpProperties,
        gain_medium=smallGainMedium,
        phi_ase=realPhiAse(crossSections),
        time_integrator=ExponentialEuler(),
        time_step_size=1e-5,
    ).on_step(seen.append)

    simulation.runSteps(2)

    assert simulation.getLastState() is seen[-1]
    assert simulation.getResults() is seen[-1]
    assert simulation.lastState is seen[-1]
    assert len(seen) == 2
    assert seen[-1].step == 2
    assert seen[-1].time == 2e-5
    assert seen[-1].betaVolume.shape == (2, 2)


def testPumpOwnsItsContributionDuration(
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    seen = []
    simulation = configuredSimulation(pumpProperties,
        gain_medium=smallGainMedium,
        phi_ase=realPhiAse(crossSections),
        time_integrator=ExponentialEuler(),
        time_step_size=1e-5,
    ).on_step(seen.append)

    simulation.step(3)

    assert np.any(seen[0].dndt_pump > 0.0)
    assert np.allclose(seen[1].dndt_pump, 0.0)
    assert np.allclose(seen[2].dndt_pump, 0.0)


def testInternalRunUsesPumpActivityWindow(
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    simulation = configuredSimulation(pumpProperties,
        gain_medium=smallGainMedium,
        phi_ase=realPhiAse(crossSections),
        time_integrator=ExponentialEuler(),
        time_step_size=1e-5,
        output_steps=(3,),
    )

    simulation.runSteps(3)

    assert pumpProperties.physical.pump_steps == 1


def testPumpRejectsNegativePumpSteps(pumpProperties):
    with pytest.raises(ValueError, match="pump_steps"):
        replace(pumpProperties.physical, pump_steps=-1)


def testSimulationDerivesStepCountFromLongestActivityWindow(
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    setup = SimpleNamespace(
        physical=replace(pumpProperties.physical, pump_steps=3),
        injector=pumpProperties.injector,
    )
    simulation = configuredSimulation(
        setup,
        gain_medium=smallGainMedium,
        phi_ase=PhiASE(spectralProperties=crossSections, openpmdBackend="adios", ase_steps=5),
        time_integrator=ExponentialEuler(),
        time_step_size=1e-5,
    )

    simulation.step()

    assert fakeCppSimulation[-1]["steps"] == 5


def testZeroActivityCountsDisablePumpAndAse(
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    setup = SimpleNamespace(
        physical=replace(pumpProperties.physical, pump_steps=0),
        injector=pumpProperties.injector,
    )
    seen = []
    simulation = configuredSimulation(
        setup,
        gain_medium=smallGainMedium,
        phi_ase=PhiASE(spectralProperties=crossSections, openpmdBackend="adios", ase_steps=0),
        time_integrator=ExponentialEuler(),
        time_step_size=1e-5,
        simulation_steps=2,
    ).on_step(seen.append)

    simulation.step()

    assert len(seen) == 2
    assert all(np.allclose(state.dndt_pump, 0.0) for state in seen)


def testDisabledActivitiesRequireExplicitSimulationSteps(
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    setup = SimpleNamespace(
        physical=replace(pumpProperties.physical, pump_steps=0),
        injector=pumpProperties.injector,
    )
    simulation = configuredSimulation(
        setup,
        gain_medium=smallGainMedium,
        phi_ase=PhiASE(spectralProperties=crossSections, openpmdBackend="adios", ase_steps=0),
        time_integrator=ExponentialEuler(),
        time_step_size=1e-5,
    )

    with pytest.raises(ValueError, match="simulation_steps is required"):
        simulation.step()


def testPumpActivityDoesNotRestartAcrossStepCalls(
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    seen = []
    simulation = configuredSimulation(
        pumpProperties,
        gain_medium=smallGainMedium,
        phi_ase=realPhiAse(crossSections),
        time_integrator=ExponentialEuler(),
        time_step_size=1e-5,
    ).on_step(seen.append)

    simulation.step(1)
    simulation.step(1)

    assert np.any(seen[0].dndt_pump > 0.0)
    assert np.allclose(seen[1].dndt_pump, 0.0)


def testStepCallbackPassesStateBeforeUserArguments(
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    seen = []

    def record(state, label, scale=1.0):
        seen.append((label, state.step, scale, state.phi_ase.shape))

    simulation = configuredSimulation(pumpProperties,
        gain_medium=smallGainMedium,
        phi_ase=realPhiAse(crossSections),
        time_integrator=ExponentialEuler(),
        time_step_size=1e-5,
    ).on_step(record, "vtk", scale=2.0)

    simulation.step(2)

    assert seen == [("vtk", 1, 2.0, (2, 2)), ("vtk", 2, 2.0, (2, 2))]


def testInitCallbacksRunBeforeCompiledTransport(
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    events = []

    def init(simulation, label, enabled=False):
        events.append(("init", label, enabled, simulation.current_step))
        simulation._testInitialized = True

    simulation = configuredSimulation(pumpProperties,
        gain_medium=smallGainMedium,
        phi_ase=realPhiAse(crossSections),
        time_integrator=ExponentialEuler(),
        time_step_size=1e-5,
    ).on_init(init, "setup", enabled=True)

    simulation.step(2)
    simulation.step(1)

    assert events == [("init", "setup", True, 0)]
    assert simulation._testInitialized is True
    assert simulation.current_step == 3
    assert simulation.current_time == 3.0000000000000004e-5


def testCompiledSimulationRejectsPythonControlCallbacks(
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    simulation = configuredSimulation(pumpProperties,
        gain_medium=smallGainMedium,
        phi_ase=realPhiAse(crossSections),
        time_integrator="explicit-euler",
        time_step_size=1e-5,
    ).before_step(lambda simulation: None)

    with pytest.raises(ValueError, match="before_step"):
        simulation.runSteps(1)


def testAutonomousFinalIsAnOutputScheduleHelper():
    assert autonomous_final(5) == (5,)
    with pytest.raises(ValueError, match="positive"):
        autonomous_final(0)


def testAutonomousOutputScheduleAndFieldsAreAppliedAtInitialization(
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    seen = []
    simulation = configuredSimulation(
        pumpProperties,
        gain_medium=smallGainMedium,
        phi_ase=realPhiAse(crossSections),
        time_integrator=ExponentialEuler(),
        time_step_size=1e-5,
        output_steps=(2, 5),
        output_fields=("beta_volume", "dndt_pump"),
    ).on_step(seen.append)

    simulation.runSteps(5)

    assert [state.step for state in seen] == [2, 5]
    assert seen[-1].betaVolume is not None
    assert seen[-1].dndtPump is not None
    assert seen[-1].phiAse is None
    assert simulation.current_step == 5


def testSynchronizedDebugExchangesSelectedControlAfterEveryNonfinalStep(
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    controlled_after = []

    def control(simulation):
        controlled_after.append(simulation.current_step)
        simulation.gainMedium.get("betaVolume").value[...] = 0.75

    simulation = configuredSimulation(
        pumpProperties,
        gain_medium=smallGainMedium,
        phi_ase=realPhiAse(crossSections),
        time_integrator=ExponentialEuler(),
        time_step_size=1e-5,
        execution_mode="synchronized-debug",
        control_fields=("beta_volume",),
    ).before_step(control)

    simulation.runSteps(3)

    assert controlled_after == [1, 2]


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"execution_mode": "autonomous", "control_fields": ("beta_volume",)}, "control_fields"),
        ({"execution_mode": "synchronized-debug", "output_steps": (1,)}, "output_steps"),
        ({"output_fields": ("point_beta",)}, "unsupported output_fields"),
        ({"execution_mode": "synchronized-debug", "control_fields": ("point_beta",)}, "unsupported control_fields"),
    ],
)
def testSimulationRunContractRejectsUnsupportedModeCombinations(
    smallGainMedium,
    crossSections,
    kwargs,
    message,
):
    with pytest.raises(ValueError, match=message):
        Simulation(
            gain_medium=smallGainMedium,
            phi_ase=realPhiAse(crossSections),
            time_integrator=ExponentialEuler(),
            time_step_size=1e-5,
            **kwargs,
        )


def testSimulationExposesOnlySnakeCaseCallbackAndRunUntilNames(
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    simulation = configuredSimulation(
        pumpProperties,
        gain_medium=smallGainMedium,
        phi_ase=realPhiAse(crossSections),
        time_integrator=ExponentialEuler(),
        time_step_size=1e-5,
    )

    assert callable(simulation.on_init)
    assert callable(simulation.on_step)
    assert callable(simulation.before_step)
    assert callable(simulation.run_until)
    assert not hasattr(simulation, "onInit")
    assert not hasattr(simulation, "onStep")
    assert not hasattr(simulation, "beforeStep")
    assert not hasattr(simulation, "runUntil")


def testCompiledSimulationRejectsExternalOpenPmdSessionOwnership(
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    simulation = configuredSimulation(pumpProperties,
        gain_medium=smallGainMedium,
        phi_ase=realPhiAse(crossSections),
        time_integrator=ExponentialEuler(),
        time_step_size=1e-5,
    )

    with pytest.raises(ValueError, match=r"owns its C\+\+ openPMD lifetime"):
        simulation.runSteps(1, openpmdSession="persistent")


def testSurfacePumpInjectorRejectsEmptyDomains():
    with pytest.raises(ValueError, match="at least one surface domain"):
        SurfacePumpInjector(())


def testPumpValidatesDedicatedRayControls(pumpProperties):
    with pytest.raises(ValueError, match="ray_count"):
        replace(pumpProperties.physical, ray_count=0)
    with pytest.raises(ValueError, match="uint32"):
        replace(pumpProperties.physical, rng_seed=2**32)


def testPhysicalPumpIsSeparateFromInjectionAndOwnsSampling(pumpProperties):
    assert pumpProperties.physical.total_power == 1.0
    assert pumpProperties.physical.spectrum.weights.tolist() == [1.0]
    assert pumpProperties.physical.angular_distribution.weights.tolist() == [1.0]
    assert pumpProperties.injector.surface_domains == (1,)
    assert pumpProperties.physical.ray_count == 256
    assert pumpProperties.physical.pump_steps == 1


def testTimeIntegrationSolverIsMandatory(
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    try:
        configuredSimulation(pumpProperties,
            gain_medium=smallGainMedium,
            phi_ase=realPhiAse(crossSections),
            time_integrator=None,
            time_step_size=1e-5,
        )
    except ValueError as exc:
        assert "time_integrator" in str(exc)
    else:
        raise AssertionError("Simulation accepted a missing timeIntegrationSolver")


def testTimeIntegrationSolversCanStepSimulation(
    fakeCppSimulation,
    pumpProperties,
    crossSections,
    smallTopology,
):
    solvers = [
        ExplicitEuler(),
        Heun(),
        Midpoint(),
        RungeKutta4(),
        FrozenPhiAseRungeKutta4(),
        ImplicitEuler(iterations=2),
        ExponentialEuler(),
    ]

    for solver in solvers:
        medium = GainMedium(topology=smallTopology).withPhysicalProperties(
            betaVolume=np.zeros((smallTopology.numberOfTriangles, smallTopology.levels - 1)),
            claddingCellTypes=np.zeros(smallTopology.numberOfTriangles, dtype=np.uint32),
            nTot=2.76e20,
            crystalTFluo=9.5e-4,
            claddingNumber=1,
            claddingAbsorption=0.0,
        )
        simulation = configuredSimulation(
            pumpProperties,
            gain_medium=medium,
            phi_ase=realPhiAse(crossSections),
            time_integrator=solver,
            time_step_size=1e-5,
        ).step()
        state = simulation.get_last_state()

        assert state.beta_volume.shape == (smallTopology.numberOfTriangles, smallTopology.levels - 1)
        assert np.all(np.isfinite(state.beta_volume))


def testPicmiStyleStepDefaultsToOneStep(
    fakeCppSimulation,
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    simulation = configuredSimulation(
        pumpProperties,
        gain_medium=smallGainMedium,
        phi_ase=realPhiAse(crossSections),
        time_integrator=ExponentialEuler(),
        time_step_size=1e-5,
        simulation_steps=2,
    )

    assert simulation.step() is simulation
    assert simulation.current_step == 2
    assert fakeCppSimulation[-1]["steps"] == 2
    assert simulation.get_last_state().step == 2


def testSimulationCanRegisterMultiplePhysicalPumps(
    smallGainMedium,
    pumpProperties,
    crossSections,
):
    simulation = configuredSimulation(
        pumpProperties,
        gain_medium=smallGainMedium,
        phi_ase=realPhiAse(crossSections),
        time_integrator=ExponentialEuler(),
        time_step_size=1e-5,
    )
    second = replace(pumpProperties.physical, total_power=2.0, name="second")

    returned = simulation.add_pump(second, injection_method=pumpProperties.injector)

    assert returned is simulation
    assert len(simulation.pump.sources) == 2
    assert [source.totalPower for source in simulation.pump.sources] == [1.0, 2.0]
