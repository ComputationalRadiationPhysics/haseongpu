#!/usr/bin/env python3
"""Launch the laserPumpCladding example with CI smoke-test defaults."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


def repoRoot() -> Path:
    return Path(os.environ.get("GITHUB_WORKSPACE", Path(__file__).resolve().parents[1]))


def launchCommand(openpmdBackend: str, outputDir: Path) -> list[str]:
    return [
        sys.executable,
        str(repoRoot() / "example" / "laserPumpCladding.py"),
        "--backend",
        "Host_Cpu_CpuSerial",
        "--openpmd-backend",
        openpmdBackend,
        "--timeSteps",
        "1",
        "--pumpSteps",
        "1",
        "--vtk-output-dir",
        str(outputDir),
        "--min-sample-range",
        "0",
        "--max-sample-range",
        "100",
        "--rng-seed",
        "5489",
    ]


def launchLaserPump(openpmdBackend: str, outputDir: Path) -> int:
    outputDir = outputDir.resolve()
    outputDir.mkdir(parents=True, exist_ok=True)
    return subprocess.run(
        launchCommand(openpmdBackend, outputDir),
        cwd=outputDir,
        check=False,
    ).returncode


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("openpmdBackend")
    parser.add_argument("outputDir", type=Path)
    args = parser.parse_args(argv)
    return launchLaserPump(args.openpmdBackend, args.outputDir)


if __name__ == "__main__":
    raise SystemExit(main())
