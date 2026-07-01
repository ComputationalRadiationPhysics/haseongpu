# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest


repoRoot = Path(__file__).resolve().parents[3]
exampleDir = repoRoot / "example"
sys.path.insert(0, str(exampleDir))
import laserPumpCladding  # noqa: E402
from pyInclude.openpmd import transport  # noqa: E402


def _vtkScalarNames(path):
    tokens = path.read_text(encoding="utf-8").split()
    return {tokens[index + 1] for index, token in enumerate(tokens) if token.upper() == "SCALARS"}


def testLaserPumpCladdingMediumUsesPrimitiveReflectivityShape():
    medium = laserPumpCladding.laserPumpCladdingMedium()

    assert medium.get("reflectivities").expectedShape == (
        medium.topology.numberOfTriangles,
        2,
    )
    assert medium.getTriangles()["reflectivities"].shape == (
        medium.topology.numberOfTriangles,
        2,
    )


@pytest.fixture
def fakeCompiledSnapshots(monkeypatch):
    calls = []

    def fake_run_simulation(simulation, *, steps, pumpSteps=None, transport=None):
        calls.append({"steps": steps, "pumpSteps": pumpSteps, "transport": transport, "phiASE": simulation.phiASE})
        point_shape = simulation.gainMedium.get("betaCells").expectedShape
        volume_shape = simulation.gainMedium.get("betaVolume").expectedShape
        states = []
        for step in range(1, steps + 1):
            pump_active = pumpSteps is None or step <= pumpSteps
            states.append(
                SimpleNamespace(
                    step=step,
                    time=step * simulation.timeStep,
                    betaCells=np.full(point_shape, 0.05 * step, dtype=np.float64),
                    betaVolume=np.full(volume_shape, 0.025 * step, dtype=np.float64),
                    phiAse=np.ones(point_shape, dtype=np.float64),
                    dndtAse=np.zeros(point_shape, dtype=np.float64),
                    dndtPump=(np.ones(point_shape, dtype=np.float64) if pump_active else np.zeros(point_shape, dtype=np.float64)),
                    aseResult=object(),
                )
            )
        return states

    monkeypatch.setattr(transport, "runSimulation", fake_run_simulation)
    return calls

def testLaserPumpCladdingExampleWritesVtkFromCompiledSnapshots(
    monkeypatch,
    tmp_path,
    smallGainMedium,
    fakeCompiledSnapshots,
):
    monkeypatch.setattr(laserPumpCladding, "laserPumpCladdingMedium", lambda **kwargs: smallGainMedium)

    state = laserPumpCladding.runExample(timeSlices=2, pumpSteps=1, vtkOutputDir=tmp_path)

    first = tmp_path / "laserPumpCladding_001.vtk"
    second = tmp_path / "laserPumpCladding_002.vtk"
    assert first.is_file()
    assert second.is_file()
    assert state.step == 2
    scalars = _vtkScalarNames(second)
    assert {"betaCells", "phiASE", "dndtAse", "dndtPump", "cladAbs"}.issubset(scalars)
    assert fakeCompiledSnapshots[-1]["pumpSteps"] == 1


def testLaserPumpCladdingExampleWiresOpenPmdBackend(
    monkeypatch,
    tmp_path,
    smallGainMedium,
    fakeCompiledSnapshots,
):
    monkeypatch.setattr(laserPumpCladding, "laserPumpCladdingMedium", lambda **kwargs: smallGainMedium)

    state = laserPumpCladding.runExample(
        timeSlices=2,
        pumpSteps=1,
        vtkOutputDir=tmp_path,
        openpmdBackend="hdf5",
    )

    assert state.step == 2
    assert fakeCompiledSnapshots[-1]["transport"] == "hdf5"
    assert fakeCompiledSnapshots[-1]["phiASE"].openpmdBackend == "hdf5"
    assert np.allclose(state.dndtPump, 0.0)
