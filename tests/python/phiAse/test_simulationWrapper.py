# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later


import argparse

from types import SimpleNamespace

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
    ])

    fromArgs = PhiASE.fromArgs(args)

    assert fromArgs.minRaysPerSample == 32
    assert fromArgs.maxRaysPerSample == 10000


def testPhiAseMpiRunUsesOpenPmdTransportMetadata(
    monkeypatch,
    smallGainMedium,
    crossSections,
    phiAseTestConfigPath,
):
    captured = {}

    def fakeRunPhiAse(phiAse, gainMedium, spectralProperties, **kwargs):
        captured["nPerNode"] = phiAse.nPerNode
        captured["numDevices"] = phiAse.numDevices
        captured["parallelMode"] = phiAse.parallelMode
        captured["gain_medium"] = gainMedium
        captured["spectral_properties"] = spectralProperties
        captured["openpmd_session"] = kwargs.get("openpmdSession")
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


def testSimulationRunStepsKeepsPersistentOpenPmdSessionUntilIntervalCompletes(monkeypatch):
    events = []
    openpmdSession = object()
    simulation = object.__new__(simulation_module.Simulation)
    simulation.phiASE = SimpleNamespace(
        openStream=lambda: events.append(("openStream",)) or openpmdSession,
        closeStream=lambda: events.append(("closeStream",)),
    )
    simulation._openpmdSession = None

    def fakeStep(self, *, openpmdSession=None):
        events.append(("step", openpmdSession))

    monkeypatch.setattr(simulation_module.Simulation, "step", fakeStep)

    simulation_module.Simulation.runSteps(
        simulation,
        2,
        openpmdSession="persistent",
    )

    assert events == [
        ("openStream",),
        ("step", openpmdSession),
        ("step", openpmdSession),
        ("closeStream",),
    ]


def testSimulationRunStepsDefaultsToPersistentOpenPmdSessionForStreamingBackend(monkeypatch):
    events = []
    openpmdSession = object()
    simulation = object.__new__(simulation_module.Simulation)
    simulation.phiASE = SimpleNamespace(
        openStream=lambda: events.append(("openStream",)) or openpmdSession,
        closeStream=lambda: events.append(("closeStream",)),
    )
    simulation._openpmdSession = None

    def fakeStep(self, *, openpmdSession=None):
        events.append(("step", openpmdSession))

    monkeypatch.setattr(
        simulation_module.transport,
        "_backend_spec",
        lambda: SimpleNamespace(streaming=True),
    )
    monkeypatch.setattr(simulation_module.Simulation, "step", fakeStep)

    simulation_module.Simulation.runSteps(simulation, 1)

    assert events == [
        ("openStream",),
        ("step", openpmdSession),
        ("closeStream",),
    ]


def testPhiAseNPerNodeLoadsFromArgsAndConfig():
    phiAse = PhiASE({"compute": {"nPerNode": 3}})
    assert phiAse.nPerNode == 3

    parser = argparse.ArgumentParser()
    PhiASE.addArguments(parser)
    args = parser.parse_args(["--n-per-node", "5"])

    fromArgs = PhiASE.fromArgs(args)
    assert fromArgs.nPerNode == 5
