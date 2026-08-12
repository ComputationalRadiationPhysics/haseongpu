#!/usr/bin/env python3
"""Launch the laserPumpCladding example with CI smoke-test defaults."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

DEFAULT_ALPAKA_BACKEND = "Host_Cpu_CpuSerial"


def repoRoot() -> Path:
    return Path(os.environ.get("GITHUB_WORKSPACE", Path(__file__).resolve().parents[1]))


def launchCommand(
    openpmdBackend: str,
    outputDir: Path,
    alpakaBackend: str = DEFAULT_ALPAKA_BACKEND,
) -> list[str]:
    command = [
        sys.executable,
        str(repoRoot() / "example" / "laserPumpCladdingApi.py"),
        "--backend",
        alpakaBackend,
        "--openpmd-backend",
        openpmdBackend,
        "--simulation-steps",
        "1",
        "--pump-steps",
        "1",
        "--ase-steps",
        "1",
        "--vtk-output-dir",
        str(outputDir),
        "--rng-seed",
        "5489",
    ]
    return command


def launchLaserPump(
    openpmdBackend: str,
    outputDir: Path,
    alpakaBackend: str = DEFAULT_ALPAKA_BACKEND,
) -> int:
    outputDir = outputDir.resolve()
    outputDir.mkdir(parents=True, exist_ok=True)
    return subprocess.run(
        launchCommand(openpmdBackend, outputDir, alpakaBackend),
        cwd=outputDir,
        check=False,
    ).returncode


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("openpmdBackend")
    parser.add_argument("outputDir", type=Path)
    parser.add_argument("--alpaka-backend", default=DEFAULT_ALPAKA_BACKEND)
    args = parser.parse_args(argv)
    return launchLaserPump(args.openpmdBackend, args.outputDir, args.alpaka_backend)


if __name__ == "__main__":
    raise SystemExit(main())
