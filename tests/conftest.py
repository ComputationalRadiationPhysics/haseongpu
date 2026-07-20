# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import os
import sys
from pathlib import Path
import copy
import importlib

repoRoot = Path(__file__).resolve().parents[1]
pythonTestPhiAseConfig = Path(__file__).parent / "data" / "cfg" / "phiAseTestConfig.yaml"
legacyPhiAseConfigFile = Path(__file__).parent / "data" / "cfg" / "legacy_config.yaml"
requiredHaseApi = (
    "GainMedium",
    "Grid",
    "MeshTopology",
    "OpenPmdBackends",
    "PhiASE",
    "PumpProperties",
    "SpectralDecomposition",
)


from openpmd_backend_matrix import openpmd_test_backends


def _resolve_import_path(entry):
    path = Path.cwd() if entry == "" else Path(entry)
    try:
        return path.resolve()
    except OSError:
        return path


def _is_under_repo(path):
    try:
        path.resolve().relative_to(repoRoot)
        return True
    except (OSError, ValueError):
        return False


sys.meta_path = [
    finder for finder in sys.meta_path if finder.__class__.__module__ != "_HASEonGPU_editable"
]


def _remove_checkout_import_paths():
    sys.path[:] = [entry for entry in sys.path if _resolve_import_path(entry) != repoRoot]


def _build_python_roots():
    candidates = [repoRoot / "build" / "python"]
    candidates.extend(
        sorted(
            repoRoot.glob("build/cp*/python"),
            key=lambda path: path.stat().st_mtime,
            reverse=True,
        )
    )
    return [path for path in candidates if path.is_dir()]


def _clear_hase_modules():
    for name in list(sys.modules):
        if name == "HASEonGPU" or name.startswith("HASEonGPU."):
            del sys.modules[name]
        elif name == "pyInclude" or name.startswith("pyInclude."):
            del sys.modules[name]


def _has_required_api(module):
    return all(hasattr(module, name) for name in requiredHaseApi)


def _import_hase_api():
    _remove_checkout_import_paths()
    _clear_hase_modules()
    try:
        module = importlib.import_module("HASEonGPU")
        if _has_required_api(module):
            return module
    except ModuleNotFoundError as err:
        if err.name != "HASEonGPU":
            raise
    _clear_hase_modules()
    sys.path[:0] = [str(repoRoot), *(str(path) for path in _build_python_roots())]
    module = importlib.import_module("HASEonGPU")
    if not _has_required_api(module):
        missing = ", ".join(name for name in requiredHaseApi if not hasattr(module, name))
        raise ImportError(f"HASEonGPU import did not expose required test API: {missing}")
    return module


_hase_api = _import_hase_api()
GainMedium = _hase_api.GainMedium
Grid = _hase_api.Grid
MeshTopology = _hase_api.MeshTopology
PhiASE = _hase_api.PhiASE
PumpProperties = _hase_api.PumpProperties
SpectralDecomposition = _hase_api.SpectralDecomposition

import numpy as np
import pytest


@pytest.fixture(scope="session", params=openpmd_test_backends())
def openPmdRuntimeBackend(request):
    return request.param


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
        reflectivities=np.zeros((smallTopology.numberOfTriangles, 2)),
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
