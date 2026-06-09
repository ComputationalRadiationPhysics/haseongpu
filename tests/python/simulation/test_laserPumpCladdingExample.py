# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import sys
from pathlib import Path

import numpy as np


repoRoot = Path(__file__).resolve().parents[3]
exampleDir = repoRoot / "example" / "python_example"
sys.path.insert(0, str(exampleDir))
import laserPumpCladding  # noqa: E402


def _vtkScalarNames(path):
    tokens = path.read_text(encoding="utf-8").split()
    return {tokens[index + 1] for index, token in enumerate(tokens) if token.upper() == "SCALARS"}


class _NoPumpSolver:
    def step(self, input, pump):
        return np.asarray(input["betaCell"], dtype=np.float64).copy()


class _FakePhiASE:
    def __init__(self, spectralProperties=None, **overrides):
        self.crossSections = spectralProperties
        self.spectralProperties = spectralProperties
        self.laserProperties = None
        self.backend = overrides.get("backend", "FakeBackend")
        self._shape = None

    def run(self, gainMedium=None, crossSections=None):
        self._shape = gainMedium.get("betaCells").expectedShape
        return self

    def getResults(self):
        class Result:
            pass

        result = Result()
        result.phiAse = np.ones(int(np.prod(self._shape)), dtype=np.float64)
        return result


def testLaserPumpCladdingExampleWritesVtkFromOnStep(monkeypatch, tmp_path, smallGainMedium):
    monkeypatch.setattr(laserPumpCladding, "BetaIntegrationGaussianSolver", lambda: _NoPumpSolver())
    monkeypatch.setattr(laserPumpCladding, "laserPumpCladdingMedium", lambda **kwargs: smallGainMedium)
    monkeypatch.setattr(
        laserPumpCladding.PhiASE,
        "fromYaml",
        classmethod(lambda cls, filename, **overrides: _FakePhiASE(**overrides)),
    )

    state = laserPumpCladding.runExample(timeSlices=2, vtkOutputDir=tmp_path)

    first = tmp_path / "laserPumpCladding_001.vtk"
    second = tmp_path / "laserPumpCladding_002.vtk"
    assert first.is_file()
    assert second.is_file()
    assert state.step == 2
    scalars = _vtkScalarNames(second)
    assert {"betaCells", "phiASE", "dndtAse", "cladAbs"}.issubset(scalars)
