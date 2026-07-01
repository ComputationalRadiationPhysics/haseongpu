# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import sys
from pathlib import Path

import numpy as np


repoRoot = Path(__file__).resolve().parents[3]
exampleDir = repoRoot / "example"
sys.path.insert(0, str(exampleDir))
import laserPumpCladding  # noqa: E402


def _vtkScalarNames(path):
    tokens = path.read_text(encoding="utf-8").split()
    return {tokens[index + 1] for index, token in enumerate(tokens) if token.upper() == "SCALARS"}


class _NoPumpSolver:
    def step(self, input, pump):
        return np.asarray(input["betaCell"], dtype=np.float64).copy()


class _ConstantPumpSolver:
    def step(self, input, pump):
        return np.asarray(input["betaCell"], dtype=np.float64) + 1.0e-6


class _FakePhiASE:
    def __init__(self, spectralProperties=None, **overrides):
        self.crossSections = spectralProperties
        self.spectralProperties = spectralProperties
        self.laserProperties = None
        self.backend = overrides.get("backend", "FakeBackend")
        self.openpmdBackend = overrides.get(
            "openpmdBackend",
            overrides.get("openpmd_backend", "adios"),
        )
        self._shape = None
        self.runInputs = []

    def run(self, gainMedium=None, crossSections=None, **kwargs):
        self._shape = gainMedium.get("betaCells").expectedShape
        self.runInputs.append(np.asarray(gainMedium.get("betaCells").value, dtype=np.float64).copy())
        return self

    def getResults(self):
        class Result:
            pass

        result = Result()
        result.phiAse = np.ones(int(np.prod(self._shape)), dtype=np.float64)
        return result

def testLaserPumpCladdingExampleWritesVtkFromOnStep(monkeypatch, tmp_path, smallGainMedium):
    monkeypatch.setattr(laserPumpCladding, "OneDimensionalZTraversal", lambda: _NoPumpSolver())
    monkeypatch.setattr(laserPumpCladding, "laserPumpCladdingMedium", lambda **kwargs: smallGainMedium)
    monkeypatch.setattr(
        laserPumpCladding.PhiASE,
        "fromYaml",
        classmethod(lambda cls, filename, **overrides: _FakePhiASE(**overrides)),
    )

    state = laserPumpCladding.runExample(timeSlices=2, pumpSteps=1, vtkOutputDir=tmp_path)

    first = tmp_path / "laserPumpCladding_001.vtk"
    second = tmp_path / "laserPumpCladding_002.vtk"
    assert first.is_file()
    assert second.is_file()
    assert state.step == 2
    scalars = _vtkScalarNames(second)
    assert {"betaCells", "phiASE", "dndtAse", "dndtPump", "cladAbs"}.issubset(scalars)


def testLaserPumpCladdingExampleSkipsAseBackendForInitialStep(monkeypatch, tmp_path, smallGainMedium):
    monkeypatch.setattr(laserPumpCladding, "OneDimensionalZTraversal", lambda: _ConstantPumpSolver())
    monkeypatch.setattr(laserPumpCladding, "laserPumpCladdingMedium", lambda **kwargs: smallGainMedium)
    phiAse = _FakePhiASE()
    monkeypatch.setattr(
        laserPumpCladding.PhiASE,
        "fromYaml",
        classmethod(lambda cls, filename, **overrides: phiAse),
    )

    state = laserPumpCladding.runExample(timeSlices=2, pumpSteps=1, vtkOutputDir=tmp_path)

    assert state.step == 2
    assert len(phiAse.runInputs) == 1
    assert np.any(phiAse.runInputs[0] > 0.0)


def testLaserPumpCladdingExamplePassesPumpStepsToSimulation(monkeypatch, tmp_path, smallGainMedium):
    monkeypatch.setattr(laserPumpCladding, "OneDimensionalZTraversal", lambda: _ConstantPumpSolver())
    monkeypatch.setattr(laserPumpCladding, "laserPumpCladdingMedium", lambda **kwargs: smallGainMedium)
    monkeypatch.setattr(
        laserPumpCladding.PhiASE,
        "fromYaml",
        classmethod(lambda cls, filename, **overrides: _FakePhiASE(**overrides)),
    )

    state = laserPumpCladding.runExample(timeSlices=2, pumpSteps=1, vtkOutputDir=tmp_path)

    assert np.allclose(state.dndtPump, 0.0)
