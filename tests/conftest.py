# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import sys
from pathlib import Path
import copy

repoRoot = Path(__file__).resolve().parents[1]
buildPythonRoot = repoRoot / "build" / "python"
pythonTestPhiAseConfig = Path(__file__).parent / "data" / "cfg" / "phiAseTestConfig.yaml"
legacyPhiAseConfigFile = Path(__file__).parent / "data" / "cfg" / "legacy_config.yaml"

sys.path.insert(0, str(repoRoot))
if buildPythonRoot.is_dir():
    sys.path.insert(0, str(buildPythonRoot))

from HASEonGPU import GainMedium, Grid, MeshTopology, PhiASE, PumpProperties, SpectralDecomposition

import numpy as np
import pytest


@pytest.fixture(scope="session")
def phiAseTestConfigPath():
    return pythonTestPhiAseConfig


@pytest.fixture(scope="session")
def phiAseTestConfig():
    import yaml

    with pythonTestPhiAseConfig.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


@pytest.fixture(scope="session")
def legacyPhiAseConfigPath():
    return legacyPhiAseConfigFile


@pytest.fixture(scope="session")
def legacyPhiAseConfig():
    import yaml

    with legacyPhiAseConfigFile.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


@pytest.fixture
def makePhiAseTestConfig(phiAseTestConfig):
    def make(**overrides):
        config = copy.deepcopy(phiAseTestConfig)
        for sectionName, sectionValues in overrides.items():
            if sectionName not in config or not isinstance(sectionValues, dict):
                config[sectionName] = sectionValues
                continue
            config[sectionName].update(sectionValues)
        return config

    return make


@pytest.fixture
def crossSections():
    return SpectralDecomposition.monochromatic(
        wavelength=940e-9,
        crossSectionAbsorption=0.01e-20,
        crossSectionEmission=0.02e-20,
    )


@pytest.fixture
def smallTopology():
    return MeshTopology.fromGrid(Grid(xExtent=1, yExtent=1, zExtent=0.5, tileSizeZ=0.25))


@pytest.fixture
def smallGainMedium(smallTopology):
    return GainMedium(topology=smallTopology).withPhysicalProperties(
        betaCells=np.zeros((4, 3)),
        claddingCellTypes=np.zeros(2, dtype=np.uint32),
        refractiveIndices=[1.8, 1.0, 1.8, 1.0],
        reflectivities=np.zeros((2, 2)),
        nTot=2.76e20,
        crystalTFluo=9.5e-4,
        claddingNumber=1,
        claddingAbsorption=0.0,
    )


@pytest.fixture
def pumpProperties(crossSections):
    return PumpProperties(
        spectralProperties=crossSections,
        intensity=16e3,
        pumpSubsteps=100,
        wavelength=940e-9,
        radiusX=1.5,
    )


@pytest.fixture
def makeFakePhiAse(monkeypatch, crossSections, phiAseTestConfigPath):
    def make(topology, **overrides):
        phiAse = PhiASE.fromYaml(phiAseTestConfigPath, spectralProperties=crossSections, **overrides)

        class Result:
            phiAse = np.zeros(topology.numberOfPoints * topology.levels)

        monkeypatch.setattr(phiAse, "run", lambda *args, **kwargs: phiAse)
        monkeypatch.setattr(phiAse, "getResults", lambda: Result())
        return phiAse

    return make
