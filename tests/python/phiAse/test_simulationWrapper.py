# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later


import argparse
from types import SimpleNamespace

import pytest

from pyInclude import PhiASE
import pyInclude.simulation as simulation_module

class DummyResult:
    phiAse = [1.0]
    mse = [0.0]
    totalRays = [4]
    dndtAse = [0.0]


def testSimulationRunUsesOpenPmdTransportAndStoresResults(
    monkeypatch,
    smallGainMedium,
    crossSections,
    phiAseTestConfigPath,
):
    captured = {}

    def fakeRunPhiAse(phiAse, gainMedium, spectralProperties, **kwargs):
        captured["phi_ase"] = phiAse
        captured["gain_medium"] = gainMedium
        captured["spectral_properties"] = spectralProperties
        captured["openpmd_session"] = kwargs.get("openpmdSession")
        return DummyResult()

    monkeypatch.setattr(simulation_module.transport, "runPhiASE", fakeRunPhiAse)

    phiAse = PhiASE.fromYaml(
        phiAseTestConfigPath,
        spectralProperties=crossSections,
        repetitions=1,
        adaptiveSteps=1,
        backend="Host_Cpu_CpuOmpBlocks",
        parallelMode="single",
        useReflections=False,
        rngSeed=1234,
    ).run(gainMedium=smallGainMedium)

    assert isinstance(phiAse.getResults(), DummyResult)
    assert captured["phi_ase"] is phiAse
    assert captured["gain_medium"] is smallGainMedium
    assert captured["spectral_properties"] is crossSections
    assert captured["openpmd_session"] is None
    assert captured["phi_ase"].minRaysPerSample == 1000
    assert captured["phi_ase"].useReflections is False
    assert captured["phi_ase"].rngSeed == 1234


def testPhiAseLoadsYamlAndArgumentOverrides(phiAseTestConfigPath):
    phiAse = PhiASE(phiAseTestConfigPath)

    assert phiAse.minRaysPerSample == 1000
    assert phiAse.maxRaysPerSample == 10000
    assert phiAse.repetitions == 1
    assert phiAse.backend == "Host_Cpu_CpuSerial"

    parser = argparse.ArgumentParser()
    PhiASE.addArguments(parser)
    args = parser.parse_args([
        "--phi-ase-config",
        str(phiAseTestConfigPath),
        "--min-rays-per-sample",
        "32",
        "--openpmd-backend",
        "adios-sst",
    ])

    fromArgs = PhiASE.fromArgs(args)

    assert fromArgs.minRaysPerSample == 32
    assert fromArgs.maxRaysPerSample == 10000
    assert fromArgs.openpmdBackend == "adios-sst"


def testPhiAseLoadsOpenPmdBackendFromConfig():
    assert PhiASE().openpmdBackend == "auto"
    assert PhiASE({"compute": {"openpmd_backend": "hdf5"}}).openpmdBackend == "hdf5"
    assert PhiASE({"compute": {"openpmdBackend": "adios-sst"}}).openpmdBackend == "adios-sst"


def testPhiAseMpiRunUsesOpenPmdTransportMetadata(
    monkeypatch,
    tmp_path,
    smallGainMedium,
    crossSections,
    phiAseTestConfigPath,
):
    captured = {}
    monkeypatch.chdir(tmp_path)

    def fakeRunPhiAse(phiAse, gainMedium, spectralProperties, **kwargs):
        captured["nPerNode"] = phiAse.nPerNode
        captured["numDevices"] = phiAse.numDevices
        captured["parallelMode"] = phiAse.parallelMode
        captured["gain_medium"] = gainMedium
        captured["spectral_properties"] = spectralProperties
        captured["openpmd_session"] = kwargs.get("openpmdSession")
        captured["command_prefix"] = kwargs.get("command_prefix")
        captured["workspace_dir"] = kwargs.get("workspace_dir")
        return DummyResult()

    monkeypatch.setattr(simulation_module.transport, "runPhiASE", fakeRunPhiAse)

    phiAse = PhiASE.fromYaml(
        phiAseTestConfigPath,
        spectralProperties=crossSections,
        parallelMode="mpi",
        numDevices=4,
        nPerNode=2,
    ).run(gainMedium=smallGainMedium)

    assert isinstance(phiAse.getResults(), DummyResult)
    assert captured["nPerNode"] == 2
    assert captured["numDevices"] == 4
    assert captured["parallelMode"] == "mpi"
    assert captured["gain_medium"] is smallGainMedium
    assert captured["spectral_properties"] is crossSections
    assert captured["openpmd_session"] is None
    assert captured["command_prefix"] == ["mpiexec", "-npernode", "2"]
    assert captured["workspace_dir"] == tmp_path / "IO" / "phiase_mpi"


def test_phiAseMpiPersistentSessionUsesConfiguredRanks(monkeypatch, tmp_path):
    captured = {}
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        simulation_module.transport,
        "openStream",
        lambda **kwargs: captured.update(kwargs) or object(),
    )

    PhiASE(parallelMode="mpi", nPerNode=3).openStream()

    assert captured["command_prefix"] == ["mpiexec", "-npernode", "3"]
    assert captured["workspace_dir"] == tmp_path / "IO" / "phiase_mpi"
    assert captured["transport"] == "auto"


def test_phiAseMpiPersistentSessionAllowsLauncherOverride(monkeypatch):
    captured = {}
    monkeypatch.setattr(
        simulation_module.transport,
        "openStream",
        lambda **kwargs: captured.update(kwargs) or object(),
    )

    PhiASE(parallelMode="mpi", nPerNode=3).openStream(command_prefix=["srun"])

    assert captured["command_prefix"] == ["srun"]


def test_phiAseMpiRejectsInvalidRanksPerNode(monkeypatch):
    monkeypatch.setattr(
        simulation_module.transport,
        "openStream",
        lambda **kwargs: pytest.fail("invalid MPI configuration must not open a transport"),
    )

    with pytest.raises(ValueError, match="nPerNode must be a positive integer"):
        PhiASE(parallelMode="mpi", nPerNode=0).openStream()


def testPhiAseRunUsesProvidedOpenPmdSession(
    monkeypatch,
    smallGainMedium,
    crossSections,
    phiAseTestConfigPath,
):
    captured = {}
    openpmdSession = object()

    def fakeRunPhiAse(phiAse, gainMedium, spectralProperties, **kwargs):
        captured["openpmdSession"] = kwargs.get("openpmdSession")
        return DummyResult()

    monkeypatch.setattr(simulation_module.transport, "runPhiASE", fakeRunPhiAse)

    PhiASE.fromYaml(
        phiAseTestConfigPath,
        spectralProperties=crossSections,
    ).run(gainMedium=smallGainMedium, openpmdSession=openpmdSession)

    assert captured["openpmdSession"] is openpmdSession


def testPhiAseRunForwardsConfiguredOpenPmdBackend(
    monkeypatch,
    smallGainMedium,
    crossSections,
):
    captured = {}

    def fakeRunPhiAse(phiAse, gainMedium, spectralProperties, **kwargs):
        captured["transport"] = kwargs.get("transport")
        captured["openpmdSession"] = kwargs.get("openpmdSession")
        return DummyResult()

    monkeypatch.setattr(simulation_module.transport, "runPhiASE", fakeRunPhiAse)

    PhiASE(
        {
            "compute": {
                "backend": "Host_Cpu_CpuSerial",
                "openpmd_backend": "hdf5",
            }
        },
        spectralProperties=crossSections,
    ).run(gainMedium=smallGainMedium)

    assert captured == {"transport": "hdf5", "openpmdSession": None}


def testPhiAsePersistentOpenPmdSessionCanBeOpenedReusedAndClosed(
    monkeypatch,
    smallGainMedium,
    crossSections,
    phiAseTestConfigPath,
):
    events = []
    openpmdSession = object()

    def fakeOpenStream(**kwargs):
        events.append(("openStream", kwargs))
        return openpmdSession

    def fakeCloseStream(session):
        events.append(("closeStream", session))

    def fakeRunPhiAse(phiAse, gainMedium, spectralProperties, **kwargs):
        events.append(("run", kwargs.get("openpmdSession")))
        return DummyResult()

    monkeypatch.setattr(simulation_module.transport, "openStream", fakeOpenStream)
    monkeypatch.setattr(simulation_module.transport, "closeStream", fakeCloseStream)
    monkeypatch.setattr(simulation_module.transport, "runPhiASE", fakeRunPhiAse)

    phiAse = PhiASE.fromYaml(
        phiAseTestConfigPath,
        spectralProperties=crossSections,
    )
    assert phiAse.openStream(transport="adios-sst") is openpmdSession
    phiAse.run(gainMedium=smallGainMedium, openpmdSession="persistent")
    phiAse.run(gainMedium=smallGainMedium, openpmdSession="persistent")
    phiAse.closeStream()

    assert events == [
        ("openStream", {"transport": "adios-sst"}),
        ("run", openpmdSession),
        ("run", openpmdSession),
        ("closeStream", openpmdSession),
    ]


def testPhiAseIntervalOpenPmdSessionUsesOneShotTransport(
    monkeypatch,
    smallGainMedium,
    crossSections,
    phiAseTestConfigPath,
):
    captured = {}

    def fakeRunPhiAse(phiAse, gainMedium, spectralProperties, **kwargs):
        captured["openpmdSession"] = kwargs.get("openpmdSession")
        return DummyResult()

    monkeypatch.setattr(simulation_module.transport, "runPhiASE", fakeRunPhiAse)

    PhiASE.fromYaml(
        phiAseTestConfigPath,
        spectralProperties=crossSections,
    ).run(gainMedium=smallGainMedium, openpmdSession="interval")

    assert captured["openpmdSession"] is None


def testSimulationRunStepsRejectsExternalOpenPmdSessionOwnership():
    simulation = object.__new__(simulation_module.Simulation)

    try:
        simulation_module.Simulation.runSteps(
            simulation,
            2,
            openpmdSession="persistent",
        )
    except ValueError as exc:
        assert "owns its C++ openPMD lifetime" in str(exc)
    else:
        raise AssertionError("compiled Simulation accepted external openPMD session ownership")


def testSimulationRunStepsPassesStreamingBackendToCompiledTransport(monkeypatch):
    captured = {}
    simulation = object.__new__(simulation_module.Simulation)
    simulation.phiASE = SimpleNamespace(openpmdBackend="adios-sst")
    simulation.pump = SimpleNamespace(getProperty=lambda name: None)
    simulation.timeStep = 1e-5
    simulation._initialized = True
    simulation._beforeStepCallbacks = []
    simulation._callbacks = []
    simulation._step = 0
    simulation._time = 0.0

    def fake_run_simulation(simulation_arg, *, steps, pumpSteps=None, transport=None):
        captured["simulation"] = simulation_arg
        captured["steps"] = steps
        captured["pumpSteps"] = pumpSteps
        captured["transport"] = transport
        return []

    monkeypatch.setattr(simulation_module.transport, "runSimulation", fake_run_simulation)

    simulation_module.Simulation.runSteps(simulation, 1)

    assert captured == {
        "simulation": simulation,
        "steps": 1,
        "pumpSteps": None,
        "transport": "adios-sst",
    }


def testPhiAseNPerNodeLoadsFromArgsAndConfig():
    phiAse = PhiASE({"compute": {"nPerNode": 3}})
    assert phiAse.nPerNode == 3

    parser = argparse.ArgumentParser()
    PhiASE.addArguments(parser)
    args = parser.parse_args(["--n-per-node", "5"])

    fromArgs = PhiASE.fromArgs(args)
    assert fromArgs.nPerNode == 5


def testPhiAseDefaultBackendUsesAvailableAlpakaBackend(monkeypatch):
    monkeypatch.setattr(simulation_module, "_preferredDefaultBackend", lambda: "Host_Cpu_CpuSerial")

    attributes = PhiASE().openPmdAttributes(numberOfSamples=1)

    assert attributes["backend"] == "Host_Cpu_CpuSerial"


def testPhiAseDefaultBackendFailureMentionsConfigure(monkeypatch):
    def fail():
        raise RuntimeError("Run `hase-configure` to generate a matching backend/openPMD setup.")

    monkeypatch.setattr(simulation_module, "_preferredDefaultBackend", fail)

    try:
        PhiASE().openPmdAttributes(numberOfSamples=1)
    except RuntimeError as exc:
        assert "hase-configure" in str(exc)
    else:
        raise AssertionError("expected default backend resolution to fail")
