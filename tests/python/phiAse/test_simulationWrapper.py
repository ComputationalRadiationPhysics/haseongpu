# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later


import argparse

import HASEonGPU_Bindings.HASEonGPU as HASEonGPU_Bindings
from HASEonGPU import PhiASE
import pyInclude.simulation as simulation_module


UNSPECIFIED_SAMPLE_RANGE = 2**32 - 1


class DummyResult:
    phiAse = [1.0]
    mse = [0.0]
    totalRays = [4]
    dndtAse = [0.0]


def testSimulationRunBuildsBindingInputsAndStoresResults(
    monkeypatch,
    smallGainMedium,
    crossSections,
    phiAseTestConfigPath,
):
    captured = {}

    def fakeCalcPhiAse(experiment, compute, hostMesh):
        captured["experiment"] = experiment
        captured["compute"] = compute
        captured["host_mesh"] = hostMesh
        return DummyResult()

    monkeypatch.setattr(HASEonGPU_Bindings, "calcPhiASE", fakeCalcPhiAse)

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
    assert captured["host_mesh"].numberOfTriangles == 2
    assert captured["host_mesh"].numberOfLevels == 3
    assert captured["compute"].minSampleRange == UNSPECIFIED_SAMPLE_RANGE
    assert captured["compute"].maxSampleRange == UNSPECIFIED_SAMPLE_RANGE
    assert captured["experiment"].minRaysPerSample == 1000
    assert captured["experiment"].useReflections is False
    assert captured["compute"].rngSeed == 1234


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


def testPhiAseMpiRunUsesLauncherWithoutBinding(
    monkeypatch,
    smallGainMedium,
    crossSections,
    phiAseTestConfigPath,
):
    captured = {}

    def failBinding(*args, **kwargs):
        raise AssertionError("binding path must not run for parallelMode='mpi'")

    def fakeMpiRun(phiAse, gainMedium, laser, laserProperties):
        captured["nPerNode"] = phiAse.nPerNode
        captured["numDevices"] = phiAse.numDevices
        captured["gain_medium"] = gainMedium
        captured["laser"] = laser
        captured["max_sigma_e"] = laserProperties.maxSigmaE
        return DummyResult()

    monkeypatch.setattr(HASEonGPU_Bindings, "calcPhiASE", failBinding)
    monkeypatch.setattr(simulation_module.mpiLauncher, "runPhiaseMPI", fakeMpiRun)

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
    assert captured["gain_medium"] is smallGainMedium
    assert captured["laser"]["l_res"] == 1
    assert captured["max_sigma_e"] == crossSections.crossSectionEmission[0]


def testPhiAseNPerNodeLoadsFromArgsAndConfig():
    phiAse = PhiASE({"compute": {"nPerNode": 3}})
    assert phiAse.nPerNode == 3

    parser = argparse.ArgumentParser()
    PhiASE.addArguments(parser)
    args = parser.parse_args(["--n-per-node", "5"])

    fromArgs = PhiASE.fromArgs(args)
    assert fromArgs.nPerNode == 5
