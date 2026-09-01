# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later


import argparse
import os
import shlex
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml
import numpy as np

from openpmd_backend_matrix import openpmd_test_backends
from alpaka_backend_matrix import alpaka_runtime_backend
from hase_units import units
from material_library import CrossSectionTable, Material

from pyInclude import Domain, GainMedium, OpticalComponent, PhiASE, VolumeTopology
import pyInclude.simulation as simulation_module


@pytest.fixture
def physicalGainMedium(crossSections):
    wavelength = float(crossSections.wavelengthsAbsorption[0])
    topology = VolumeTopology.fromTetrahedra(
        np.eye(4, 3),
        [[0, 1, 2, 3]],
    )
    material = Material(
        materialName="test gain material",
        temperature=293.15 * units.K,
        refractiveIndex=1.8,
        fluorescenceLifetime=9.5e-4 * units.s,
        crossSections=CrossSectionTable.monochromatic(
            wavelength=wavelength * units.m,
            absorption=crossSections.absorptionAt(wavelength) * units.cm**2,
            emission=crossSections.emissionAt(wavelength) * units.cm**2,
        ),
        active=True,
        activeIonDensity=2.76e20 / units.cm**3,
    )
    return GainMedium(
        [OpticalComponent(domain=Domain.fromTopology(topology), material=material)]
    )

class DummyResult:
    phiAse = [1.0]
    standardError = [0.0]
    relativeStandardError = [0.0]
    totalRays = [4]
    dndtAse = [0.0]


def testCiAlpakaRuntimeBackendSelectionUsesCudaWhenRequested(monkeypatch):
    backends = ["Host_Cpu_CpuOmpBlocks", "Nvidia_Gpu_CudaRt"]
    monkeypatch.setenv("HASE_TEST_ALPAKA_BACKEND", "cuda")

    assert alpaka_runtime_backend(backends) == "Nvidia_Gpu_CudaRt"


def testCiAlpakaRuntimeBackendSelectionRejectsMissingCuda(monkeypatch):
    monkeypatch.setenv("HASE_TEST_ALPAKA_BACKEND", "cuda")

    with pytest.raises(RuntimeError, match="exactly one CUDA GPU Alpaka backend"):
        alpaka_runtime_backend(["Host_Cpu_CpuOmpBlocks"])


def testCiAlpakaRuntimeBackendSelectionAcceptsExactEnvironmentName(monkeypatch):
    requested = "Nvidia_Gpu_CudaRt"
    monkeypatch.setenv("HASE_TEST_ALPAKA_BACKEND", requested)

    assert alpaka_runtime_backend(["Host_Cpu_CpuSerial", requested]) == requested


def testSimulationRunUsesOpenPmdTransportAndStoresResults(
    monkeypatch,
    physicalGainMedium,
    crossSections,
    phiAseTestConfigPath,
):
    captured = {}

    def fakeRunPhiAse(request, **kwargs):
        captured["phi_ase"] = request.phiASE
        captured["gain_medium"] = request.gainMedium
        captured["openpmd_session"] = kwargs.get("openpmdSession")
        return DummyResult()

    monkeypatch.setattr(simulation_module.transport, "runPhiASE", fakeRunPhiAse)

    phiAse = PhiASE.fromYaml(
        phiAseTestConfigPath,
        repetitions=1,
        adaptiveSteps=1,
        parallelMode="single",
        useReflections=False,
        rngSeed=1234,
    ).run(gainMedium=physicalGainMedium)

    assert isinstance(phiAse.getResults(), DummyResult)
    assert captured["phi_ase"] is phiAse
    assert captured["gain_medium"] is physicalGainMedium
    transportedMaterial = captured["gain_medium"].components[0].material
    assert transportedMaterial.crossSections.wavelengths.toValue(units.m)[0] == pytest.approx(
        crossSections.wavelengthsAbsorption[0]
    )
    assert transportedMaterial.crossSections.absorption.toValue(units.cm**2)[0] == pytest.approx(
        crossSections.crossSectionAbsorption[0]
    )
    assert transportedMaterial.crossSections.emission.toValue(units.cm**2)[0] == pytest.approx(
        crossSections.crossSectionEmission[0]
    )
    assert captured["openpmd_session"] is None
    assert captured["phi_ase"].minRays == 1000
    assert captured["phi_ase"].useReflections is False
    assert captured["phi_ase"].rngSeed == 1234


def testPhiAseRejectsUnavailableAlpakaBackendBeforeTransport(
    monkeypatch,
    physicalGainMedium,
    crossSections,
):
    monkeypatch.setattr(
        simulation_module.AlpakaBackends,
        "all",
        lambda: ["Host_Cpu_CpuSerial", "Host_NumaCpu_CpuSerial"],
    )
    monkeypatch.setattr(
        simulation_module.transport,
        "_ensure_backend_available",
        lambda backend: pytest.fail("openPMD preflight must follow compute preflight"),
    )
    monkeypatch.setattr(
        simulation_module.transport,
        "runPhiASE",
        lambda *args, **kwargs: pytest.fail("an invalid compute backend must not start transport"),
    )

    phi_ase = PhiASE(
        backend="Host_Cpu_CpuOmpBlocks",
    )
    with pytest.raises(RuntimeError, match="Host_Cpu_CpuOmpBlocks") as error:
        phi_ase.run(gainMedium=physicalGainMedium)

    assert "Host_Cpu_CpuSerial" in str(error.value)
    assert "Host_NumaCpu_CpuSerial" in str(error.value)


def testPhiAseRejectsUnavailableOpenPmdBackendBeforeTransport(
    monkeypatch,
    physicalGainMedium,
    crossSections,
):
    monkeypatch.setattr(simulation_module.AlpakaBackends, "all", lambda: ["Host_Cpu_CpuSerial"])

    def reject_openpmd(backend):
        raise RuntimeError(
            f"openPMD backend '{backend}' is unavailable; available backends: adios, adios-sst"
        )

    monkeypatch.setattr(simulation_module.transport, "_ensure_backend_available", reject_openpmd)
    monkeypatch.setattr(
        simulation_module.transport,
        "runPhiASE",
        lambda *args, **kwargs: pytest.fail("an invalid openPMD backend must not start transport"),
    )

    phi_ase = PhiASE(
        backend="Host_Cpu_CpuSerial",
        openpmdBackend="hdf5",
    )
    with pytest.raises(RuntimeError, match="available backends: adios, adios-sst"):
        phi_ase.run(gainMedium=physicalGainMedium)


def testPhiAseMpiPreflightDoesNotUseLauncherLocalDeviceVisibility(monkeypatch):
    monkeypatch.setattr(
        simulation_module.AlpakaBackends,
        "all",
        lambda: pytest.fail("MPI compute availability is rank-local"),
    )
    checked = []
    monkeypatch.setattr(
        simulation_module.transport,
        "_ensure_backend_available",
        lambda backend: checked.append(backend),
    )

    simulation_module._validate_launch_backends(
        PhiASE(
            backend="Cuda_NvidiaGpu_GpuCuda",
            openpmdBackend="adios",
            parallelMode="mpi",
        )
    )

    assert checked == ["adios"]


def testPhiAseLoadsYamlAndArgumentOverrides(phiAseTestConfigPath):
    phiAse = PhiASE.fromYaml(phiAseTestConfigPath)

    assert phiAse.minRays == 1000
    assert phiAse.maxRays == 10000
    assert phiAse.repetitions == 1
    assert phiAse.backend == "Host_Cpu_CpuSerial"
    assert phiAse.reflectionMode == "direct"
    assert phiAse.surfaceReservoirSize == 64
    assert phiAse.srmPositionMode == "exact"
    assert phiAse.trackRayVisits is False

    parser = argparse.ArgumentParser()
    PhiASE.addArguments(parser)
    args = parser.parse_args([
        "--phi-ase-config",
        str(phiAseTestConfigPath),
        "--min-rays",
        "32",
        "--track-ray-visits",
        "--openpmd-backend",
        "adios-sst",
        "--use-reflections",
        "--reflection-mode",
        "srm",
        "--surface-reservoir-size",
        "256",
        "--srm-position-mode",
        "centroid",
        "--monochromatic",
        "--write-vtk",
        "--devices",
        "2",
        "4",
        "--min-sample-range",
        "3",
        "--max-sample-range",
        "11",
        "--ase-steps",
        "7",
    ])

    fromArgs = PhiASE.fromArgs(args)

    assert fromArgs.minRays == 32
    assert fromArgs.trackRayVisits is True
    assert fromArgs.maxRays == 10000
    assert fromArgs.openpmdBackend == "adios-sst"
    assert fromArgs.useReflections is True
    assert fromArgs.reflectionMode == "srm"
    assert fromArgs.surfaceReservoirSize == 256
    assert fromArgs.srmPositionMode == "centroid"
    assert fromArgs.monochromatic is True
    assert fromArgs.writeVtk is True
    assert fromArgs.devices == [2, 4]
    assert fromArgs.minSampleRange == 3
    assert fromArgs.maxSampleRange == 11
    assert fromArgs.ase_steps == 7


@pytest.mark.parametrize(
    ("keyword", "value", "message"),
    (
        ("reflectionMode", "histogram", "reflectionMode"),
        ("surfaceReservoirSize", 0, "surfaceReservoirSize"),
        ("surfaceReservoirSize", 2.5, "surfaceReservoirSize"),
        ("srmPositionMode", "vertex", "srmPositionMode"),
    ),
)
def testPhiAseRejectsInvalidSrmConfiguration(keyword, value, message):
    with pytest.raises((TypeError, ValueError), match=message):
        PhiASE(**{keyword: value})


def testPhiAseLoadsCentroidSrmFromYaml(tmp_path):
    path = tmp_path / "centroid-srm.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "schema_version": 3,
                "simulation": {
                    "phi_ase": {
                        "use_reflections": True,
                        "reflection_mode": "srm",
                        "surface_reservoir_size": 512,
                        "srm_position_mode": "centroid",
                    }
                },
            }
        ),
        encoding="utf-8",
    )

    phiAse = PhiASE.fromYaml(path)

    assert phiAse.useReflections is True
    assert phiAse.reflectionMode == "srm"
    assert phiAse.surfaceReservoirSize == 512
    assert phiAse.srmPositionMode == "centroid"

def testPhiAseDefaultBackendSerializesAvailableAlpakaBackend():
    phiAse = PhiASE()
    backend = phiAse.openPmdAttributes(numberOfSamples=1)["backend"]

    assert backend != "gpu"
    assert backend in simulation_module.AlpakaBackends.all()


def testPhiAseSerializesAdaptiveRangeWithoutAnImplicitFixedRayCount():
    phiAse = PhiASE(minRays=100, maxRays=1600, adaptiveSteps=4)

    attributes = phiAse.openPmdAttributes(numberOfSamples=1)

    assert attributes["minRays"] == 100
    assert attributes["maxRays"] == 1600
    assert attributes["adaptiveSteps"] == 4
    assert attributes["forwardRayCount"] == 0
    assert attributes["trackRayVisits"] is False

    fixed = PhiASE(minRays=100, maxRays=1600, forwardRayCount=250, trackRayVisits=True)
    assert fixed.openPmdAttributes(numberOfSamples=1)["forwardRayCount"] == 250
    assert fixed.openPmdAttributes(numberOfSamples=1)["trackRayVisits"] is True

def testPhiAseDefaultBackendSerializesAvailableAlpakaBackend():
    phiAse = PhiASE()
    backend = phiAse.openPmdAttributes(numberOfSamples=1)["backend"]

    assert backend != "gpu"
    assert backend in simulation_module.AlpakaBackends.all()


def testPhiAseLoadsOpenPmdBackendFromConfig():
    assert PhiASE().openpmdBackend == "auto"
    assert PhiASE(openpmdBackend="hdf5").openpmdBackend == "hdf5"
    assert PhiASE(openpmdBackend="adios-sst").openpmdBackend == "adios-sst"


def testPhiAseRejectsRetiredMseThreshold(tmp_path):
    path = tmp_path / "retired-mse.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "schema_version": 3,
                "simulation": {"phi_ase": {"mse_threshold": 0.1}},
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="mse_threshold"):
        PhiASE.fromYaml(path)


def testPhiAseMpiRunUsesOpenPmdTransportMetadata(
    monkeypatch,
    tmp_path,
    physicalGainMedium,
    crossSections,
    phiAseTestConfigPath,
):
    captured = {}
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("HASE_MPIEXEC_EXTRA_ARGS", raising=False)

    def fakeRunPhiAse(request, **kwargs):
        captured["nPerNode"] = request.phiASE.nPerNode
        captured["numDevices"] = request.phiASE.numDevices
        captured["parallelMode"] = request.phiASE.parallelMode
        captured["gain_medium"] = request.gainMedium
        captured["openpmd_session"] = kwargs.get("openpmdSession")
        captured["command_prefix"] = kwargs.get("command_prefix")
        captured["workspace_dir"] = kwargs.get("workspace_dir")
        return DummyResult()

    monkeypatch.setattr(simulation_module.transport, "runPhiASE", fakeRunPhiAse)

    phiAse = PhiASE.fromYaml(
        phiAseTestConfigPath,
        parallelMode="mpi",
        numDevices=4,
        nPerNode=2,
    ).run(gainMedium=physicalGainMedium)

    assert isinstance(phiAse.getResults(), DummyResult)
    assert captured["nPerNode"] == 2
    assert captured["numDevices"] == 4
    assert captured["parallelMode"] == "mpi"
    assert captured["gain_medium"] is physicalGainMedium
    assert captured["gain_medium"].components[0].material.crossSections is not None
    assert captured["openpmd_session"] is None
    assert captured["command_prefix"] == ["mpiexec", "-npernode", "2"]
    assert captured["workspace_dir"] == tmp_path / "IO" / "phiase_mpi"


def test_phiAseMpiPersistentSessionUsesConfiguredRanks(monkeypatch, tmp_path):
    captured = {}
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("HASE_MPIEXEC_EXTRA_ARGS", raising=False)
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
    physicalGainMedium,
    crossSections,
    phiAseTestConfigPath,
):
    captured = {}
    openpmdSession = object()

    def fakeRunPhiAse(request, **kwargs):
        captured["openpmdSession"] = kwargs.get("openpmdSession")
        return DummyResult()

    monkeypatch.setattr(simulation_module.transport, "runPhiASE", fakeRunPhiAse)

    PhiASE.fromYaml(
        phiAseTestConfigPath,
    ).run(gainMedium=physicalGainMedium, openpmdSession=openpmdSession)

    assert captured["openpmdSession"] is openpmdSession


def testPhiAseRunForwardsConfiguredOpenPmdBackend(
    monkeypatch,
    physicalGainMedium,
    crossSections,
):
    captured = {}
    preflight = []

    def fakeRunPhiAse(request, **kwargs):
        captured["transport"] = kwargs.get("transport")
        captured["openpmdSession"] = kwargs.get("openpmdSession")
        return DummyResult()

    monkeypatch.setattr(simulation_module.transport, "runPhiASE", fakeRunPhiAse)
    monkeypatch.setattr(
        simulation_module.transport,
        "_ensure_backend_available",
        lambda backend: preflight.append(backend),
    )

    PhiASE(
        backend="Host_Cpu_CpuSerial",
        openpmdBackend="hdf5",
    ).run(gainMedium=physicalGainMedium)

    assert captured == {"transport": "hdf5", "openpmdSession": None}
    assert preflight == ["hdf5"]


def testPhiAsePersistentOpenPmdSessionCanBeOpenedReusedAndClosed(
    monkeypatch,
    physicalGainMedium,
    crossSections,
    phiAseTestConfigPath,
):
    if "adios-sst" not in openpmd_test_backends():
        pytest.skip("persistent openPMD session test requires the adios-sst backend")

    events = []
    openpmdSession = object()

    def fakeOpenStream(**kwargs):
        events.append(("openStream", kwargs))
        return openpmdSession

    def fakeCloseStream(session):
        events.append(("closeStream", session))

    def fakeRunPhiAse(request, **kwargs):
        events.append(("run", kwargs.get("openpmdSession")))
        return DummyResult()

    monkeypatch.setattr(simulation_module.transport, "openStream", fakeOpenStream)
    monkeypatch.setattr(simulation_module.transport, "closeStream", fakeCloseStream)
    monkeypatch.setattr(simulation_module.transport, "runPhiASE", fakeRunPhiAse)

    phiAse = PhiASE.fromYaml(
        phiAseTestConfigPath,
    )
    assert phiAse.openStream(transport="adios-sst") is openpmdSession
    phiAse.run(gainMedium=physicalGainMedium, openpmdSession="persistent")
    phiAse.run(gainMedium=physicalGainMedium, openpmdSession="persistent")
    phiAse.closeStream()

    expected_open_kwargs = {"transport": "adios-sst"}
    if phiAse.parallelMode == "mpi":
        expected_open_kwargs.update(
            {
                "command_prefix": [
                    "mpiexec",
                    *shlex.split(os.environ.get("HASE_MPIEXEC_EXTRA_ARGS", "")),
                    "-npernode",
                    str(phiAse.nPerNode),
                ],
                "workspace_dir": Path.cwd() / "IO" / "phiase_mpi",
            }
        )

    assert events == [
        ("openStream", expected_open_kwargs),
        ("run", openpmdSession),
        ("run", openpmdSession),
        ("closeStream", openpmdSession),
    ]


def testPhiAseIntervalOpenPmdSessionUsesOneShotTransport(
    monkeypatch,
    physicalGainMedium,
    crossSections,
    phiAseTestConfigPath,
):
    captured = {}

    def fakeRunPhiAse(request, **kwargs):
        captured["openpmdSession"] = kwargs.get("openpmdSession")
        return DummyResult()

    monkeypatch.setattr(simulation_module.transport, "runPhiASE", fakeRunPhiAse)

    PhiASE.fromYaml(
        phiAseTestConfigPath,
    ).run(gainMedium=physicalGainMedium, openpmdSession="interval")

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
    if "adios-sst" not in openpmd_test_backends():
        pytest.skip("streaming simulation transport test requires the adios-sst backend")

    captured = {}
    simulation = object.__new__(simulation_module.Simulation)
    simulation.phiASE = SimpleNamespace(
        openpmdBackend="adios-sst",
        _transportLaunchOptions=lambda: {},
    )
    simulation.pump = SimpleNamespace(getProperty=lambda name: None)
    simulation.timeStep = 1e-5
    simulation._initialized = True
    simulation._before_step_callbacks = []
    simulation._step_callbacks = []
    simulation._step = 0
    simulation._time = 0.0
    simulation.reportTimings = False

    def fake_run_simulation(simulation_arg, *, steps, transport=None, on_state=None):
        captured["simulation"] = simulation_arg
        captured["steps"] = steps
        captured["transport"] = transport
        return []

    monkeypatch.setattr(simulation_module.transport, "runSimulation", fake_run_simulation)

    simulation_module.Simulation.runSteps(simulation, 1)

    assert captured == {
        "simulation": simulation,
        "steps": 1,
        "transport": "adios-sst",
    }


def testPhiAseNPerNodeLoadsFromArgsAndConfig():
    phiAse = PhiASE(nPerNode=3)
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
