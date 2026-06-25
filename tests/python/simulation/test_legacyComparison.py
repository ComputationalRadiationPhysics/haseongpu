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
from pyInclude import AlpakaBackends, GainMedium

pytest.importorskip("scipy")


repoRoot = Path(__file__).resolve().parents[3]
legacyScript = repoRoot / "example" / "python_example" / "legacy" / "laserPumpCladdingExample.py"
legacyMaterial = repoRoot / "example" / "python_example" / "legacy" / "pt.mat"
newExampleDir = repoRoot / "example" / "python_example"
sys.path.insert(0, str(newExampleDir))
from laserPumpCladdingOldPump import runExample  # noqa: E402


alpakaBackends = AlpakaBackends.all()
REGRESSION_RNG_SEED = 5489


def _runLegacyScript(workdir, phiAseConfig, backend, minSample_i, maxSample_i, rngSeed=None):
    experiment = phiAseConfig["experiment"]
    compute = phiAseConfig["compute"]
    command = [
        sys.executable,
        str(legacyScript),
        "--material",
        str(legacyMaterial),
        "--backend",
        str(backend),
        "--number_of_gpus",
        str(compute["numDevices"]),
        "--parallel_mode",
        str(compute["parallelMode"]),
        "--min-rays",
        str(experiment["minRaysPerSample"]),
        "--max-rays",
        str(experiment["maxRaysPerSample"]),
        "--reflection",
        str(int(bool(experiment["useReflections"]))),
        "--repetitions",
        str(compute["repetitions"]),
        "--adaptive-steps",
        str(compute["adaptiveSteps"]),
        "--mse-threshold",
        str(experiment["mseThreshold"]),
        "--timeslice",
        "50",
        "--timeslice-total",
        "2",
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
def testTimeSteppedSimulationMatchesLaserPumpCladdingScript(
    backend,
    tmp_path,
    legacyPhiAseConfig,
    legacyPhiAseConfigPath,
):
    minSample_i=0
    maxSample_i=100
    timeSlices=2
    _runLegacyScript(
        tmp_path,
        legacyPhiAseConfig,
        backend,
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
        rngSeed=REGRESSION_RNG_SEED,
    )
    medium=GainMedium.fromVtk(tmp_path / f'dndt_ASE_{timeSlices}.vtk')
    legacyDndtAse=medium.dntdAse
    modernPhiAse = np.asarray(modernState.phiAse).reshape(-1, order="F")
    modernDndtAse = np.asarray(modernState.dndtAse).reshape(-1, order="F")

    assert np.isfinite(modernState.betaCells).all()
    assert np.isfinite(modernState.betaVolume).all()
    assert np.isfinite(modernPhiAse).all()
    assert len(modernPhiAse) == len(legacyDndtAse)

    assert np.allclose(modernDndtAse, legacyDndtAse, rtol=1e-5, atol=1e-4)
    assert np.all(modernState.betaCells >= 0.0)
    assert np.all(modernState.betaCells <= 1.0)
