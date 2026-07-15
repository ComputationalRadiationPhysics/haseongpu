#!/usr/bin/env python3
"""Benchmark fixed-ray PhiASE runs for a Tet4 VTK input."""

from __future__ import annotations

import argparse
import csv
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_DIR = REPO_ROOT / "example"
sys.path[:0] = [str(REPO_ROOT), str(EXAMPLE_DIR)]

from HASEonGPU import GainMedium, PhiASE  # noqa: E402
from laserPumpCladding import laserPumpCladdingSpectralProperties  # noqa: E402


DEFAULT_BACKENDS = (
    "Cuda_NvidiaGpu_GpuCuda",
    "Host_Cpu_CpuOmpBlocks",
    "Host_Cpu_CpuSerial",
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vtk", type=Path, help="Tet4 VTK input containing the gain-medium fields")
    parser.add_argument(
        "--rays",
        type=int,
        nargs="+",
        default=(1_000_000, 10_000_000, 50_000_000),
        help="Fixed global ray counts (default: 1000000 10000000 50000000)",
    )
    parser.add_argument("--backends", nargs="+", default=DEFAULT_BACKENDS)
    parser.add_argument(
        "--config",
        type=Path,
        default=REPO_ROOT / "config" / "hase-phiase.yaml",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "benchmarks" / "vtk_backend_benchmark.csv",
    )
    parser.add_argument("--spectral-resolution", type=int, default=1000)
    parser.add_argument("--rng-seed", type=int, default=12345)
    return parser.parse_args()


def main():
    args = parse_args()
    vtk_path = args.vtk.resolve()
    output_path = args.output.resolve()
    config_path = args.config.resolve()
    for path in (vtk_path, config_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    if any(rays <= 0 for rays in args.rays):
        raise ValueError("all ray counts must be positive")

    medium = GainMedium.fromVtk(vtk_path)
    spectral = laserPumpCladdingSpectralProperties(args.spectral_resolution)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fields = (
        "timestamp_utc",
        "host",
        "vtk",
        "backend",
        "rays",
        "adaptive_steps",
        "rng_seed",
        "elapsed_seconds",
        "seconds_per_million_rays",
        "million_rays_per_second",
    )
    rows = []
    for backend in args.backends:
        for rays in args.rays:
            phi_ase = PhiASE.fromYaml(
                config_path,
                spectralProperties=spectral,
                backend=backend,
                minRays=rays,
                maxRays=rays,
                adaptiveSteps=0,
                rngSeed=args.rng_seed,
            )
            started = time.perf_counter()
            phi_ase.run(gainMedium=medium, crossSections=spectral)
            elapsed = time.perf_counter() - started
            normalized = elapsed * 1_000_000 / rays
            row = {
                "timestamp_utc": datetime.now(timezone.utc).isoformat(),
                "host": platform.node(),
                "vtk": str(vtk_path),
                "backend": backend,
                "rays": rays,
                "adaptive_steps": 0,
                "rng_seed": args.rng_seed,
                "elapsed_seconds": f"{elapsed:.9f}",
                "seconds_per_million_rays": f"{normalized:.9f}",
                "million_rays_per_second": f"{1.0 / normalized:.9f}",
            }
            rows.append(row)
            print(
                f"{backend:30s} rays={rays:>9d} elapsed={elapsed:9.3f}s "
                f"normalized={normalized:9.3f}s/Mray",
                flush=True,
            )

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
