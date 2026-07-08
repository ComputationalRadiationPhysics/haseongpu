# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later


import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from pyInclude import AlpakaBackends
from pyInclude.geometry.vtk import _parseVtk

pytest.importorskip("scipy")


repoRoot = Path(__file__).resolve().parents[3]
exampleScript = repoRoot / "example" / "laserPumpCladding.py"
exampleDir = repoRoot / "example"
sys.path.insert(0, str(exampleDir))
from laserPumpCladding import runExample  # noqa: E402


alpakaBackends = AlpakaBackends.all()
REGRESSION_RNG_SEED = 5489


def _runExampleScript(workdir, phiAseConfigPath, backend, openPmdBackend, timeSlices, minSample_i, maxSample_i, rngSeed=None):
    command = [
        sys.executable,
        str(exampleScript),
        "--phi-ase-config",
        str(phiAseConfigPath),
        "--backend",
        str(backend),
        "--openpmd-backend",
        str(openPmdBackend),
        "--timeSteps",
        str(timeSlices),
        "--vtk-output-dir",
        str(workdir),
        "--min-sample-range",
        str(minSample_i),
        "--max-sample-range",
        str(maxSample_i),
    ]
    if rngSeed is not None:
        command.extend(["--rng-seed", str(rngSeed)])
    completed = subprocess.run(
        command,
        cwd=workdir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if completed.returncode != 0:
        pytest.fail(completed.stdout)


@pytest.mark.integration
@pytest.mark.parametrize("backend", alpakaBackends)
def test_timeSteppedSimulationMatchesLaserPumpCladdingScript(
    backend,
    tmp_path,
    legacyPhiAseConfigPath,
    openPmdRuntimeBackend,
):
    minSample_i=0
    maxSample_i=100
    timeSlices=2
    _runExampleScript(
        tmp_path,
        legacyPhiAseConfigPath,
        backend,
        openPmdRuntimeBackend,
        timeSlices,
        minSample_i,
        maxSample_i,
        rngSeed=REGRESSION_RNG_SEED,
    )
    modernState = runExample(
        phiAseConfigPath=legacyPhiAseConfigPath,
        timeSlices=timeSlices,
        minSampleRange=minSample_i,
        maxSampleRange=maxSample_i,
        backend=backend,
        openpmdBackend=openPmdRuntimeBackend,
        rngSeed=REGRESSION_RNG_SEED,
    )
    _, _, _, pointData, _, _ = _parseVtk(tmp_path / f"laserPumpCladding_{timeSlices:03d}.vtk")
    scriptDndtAse = np.asarray(pointData["dndtAse"])
    modernPhiAse = np.asarray(modernState.phiAse).reshape(-1, order="F")
    modernDndtAse = np.asarray(modernState.dndtAse).reshape(-1, order="F")

    assert np.isfinite(modernState.betaCells).all()
    assert np.isfinite(modernState.betaVolume).all()
    assert np.isfinite(modernPhiAse).all()
    assert len(modernPhiAse) == len(scriptDndtAse)

    assert np.allclose(modernDndtAse, scriptDndtAse, rtol=1e-5, atol=1e-4)
    assert np.all(modernState.betaCells >= 0.0)
    assert np.all(modernState.betaCells <= 1.0)
