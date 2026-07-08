# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import csv
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest


repoRoot = Path(__file__).resolve().parents[3]
exampleScript = repoRoot / "example" / "laserPumpCladding.py"
plotSsgScript = repoRoot / "scripts" / "plot_ssg.py"
juliaReferencePath = repoRoot / "tests" / "data" / "julia1D.csv"
TIME_STEP_S = 2.0e-5


def _runCommand(command, *, cwd, env=None):
    completed = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if completed.returncode != 0:
        pytest.fail(completed.stdout)


def _readColumns(path, columns):
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {
        column: np.asarray([float(row[column]) for row in rows], dtype=np.float64)
        for column in columns
    }


@pytest.mark.integration
def test_julia1DMatchesDisabledAse(tmp_path):
    vtkOutputDir = tmp_path / "vtk"
    plotPrefix = tmp_path / "ssg_z_origin"
    vtkOutputDir.mkdir()

    _runCommand(
        [
            sys.executable,
            str(exampleScript),
            "--disable-ase",
            "--timeSteps",
            "100",
            "--pumpSteps",
            "50",
            "--vtk-output-dir",
            str(vtkOutputDir),
        ],
        cwd=repoRoot,
    )

    env = os.environ.copy()
    env.setdefault("MPLBACKEND", "Agg")
    _runCommand(
        [
            sys.executable,
            str(plotSsgScript),
            "--input-dir",
            str(vtkOutputDir),
            "--direction",
            "z",
            "--x",
            "0",
            "--y",
            "0",
            "--field",
            "localGain",
            "--output-prefix",
            str(plotPrefix),
        ],
        cwd=repoRoot,
        env=env,
    )

    reference = _readColumns(juliaReferencePath, ("step", "time_s", "SSG"))
    plotted = _readColumns(
        plotPrefix.with_suffix(".csv"),
        ("timestep", "net_gain_factor"),
    )

    assert plotPrefix.with_suffix(".png").is_file()
    assert reference["step"][0] == 0.0
    assert reference["time_s"][0] == 0.0

    reference_slice = slice(1, None)
    plotted_time_s = plotted["timestep"] * TIME_STEP_S
    np.testing.assert_array_equal(plotted["timestep"], reference["step"][reference_slice])
    np.testing.assert_allclose(
        plotted_time_s,
        reference["time_s"][reference_slice],
        rtol=1.0e-4,
        atol=0.0,
    )
    np.testing.assert_allclose(
        plotted["net_gain_factor"],
        reference["SSG"][reference_slice],
        rtol=1.0e-4,
        atol=0.0,
    )
