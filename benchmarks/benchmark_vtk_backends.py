#!/usr/bin/env python3
"""Benchmark fixed-ray PhiASE runs for a Tet4 VTK input."""

from __future__ import annotations

import argparse
import csv
import platform
import statistics
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
        "--openpmd-backend",
        help="Override the openPMD storage backend selected by the YAML configuration.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "benchmarks" / "vtk_backend_benchmark.csv",
    )
    parser.add_argument("--spectral-resolution", type=int, default=1000)
    parser.add_argument("--rng-seed", type=int, default=12345)
    parser.add_argument("--warmups", type=int, default=0)
    parser.add_argument("--repetitions", type=int, default=1)
    parser.add_argument(
        "--reflection-mode",
        choices=("direct", "reflection", "both"),
        default="direct",
        help="Benchmark direct propagation, reflection-enabled propagation, or both.",
    )
    parser.add_argument("--reflection-max-iterations", type=int, default=2)
    parser.add_argument("--reflection-tolerance", type=float, default=0.0)
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
    if args.warmups < 0:
        raise ValueError("warmups must be non-negative")
    if args.repetitions <= 0:
        raise ValueError("repetitions must be positive")
    if args.reflection_max_iterations <= 0:
        raise ValueError("reflection-max-iterations must be positive")
    if args.reflection_tolerance < 0.0:
        raise ValueError("reflection-tolerance must be non-negative")

    medium = GainMedium.fromVtk(vtk_path)
    spectral = laserPumpCladdingSpectralProperties(args.spectral_resolution)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fields = (
        "timestamp_utc",
        "host",
        "vtk",
        "backend",
        "openpmd_backend",
        "rays",
        "mode",
        "repetition",
        "adaptive_steps",
        "rng_seed",
        "reflection_max_iterations",
        "reflection_tolerance",
        "elapsed_seconds",
        "seconds_per_million_rays",
        "million_rays_per_second",
    )
    rows = []
    modes = ("direct", "reflection") if args.reflection_mode == "both" else (args.reflection_mode,)
    for backend in args.backends:
        for rays in args.rays:
            for mode in modes:
                use_reflections = mode == "reflection"

                def run_once():
                    phi_ase = PhiASE.fromYaml(
                        config_path,
                        spectralProperties=spectral,
                        backend=backend,
                        **({"openpmdBackend": args.openpmd_backend} if args.openpmd_backend else {}),
                        minRays=rays,
                        maxRays=rays,
                        forwardRayCount=rays,
                        repetitions=1,
                        adaptiveSteps=0,
                        useReflections=use_reflections,
                        reflectionMaxIterations=args.reflection_max_iterations,
                        reflectionTolerance=args.reflection_tolerance,
                        rngSeed=args.rng_seed,
                    )
                    started = time.perf_counter()
                    phi_ase.run(gainMedium=medium, crossSections=spectral)
                    return time.perf_counter() - started, phi_ase.openpmdBackend

                for warmup in range(args.warmups):
                    elapsed, _openpmd_backend = run_once()
                    print(
                        f"warmup {warmup + 1}/{args.warmups}: {backend:30s} mode={mode:10s} "
                        f"rays={rays:>9d} elapsed={elapsed:9.3f}s",
                        flush=True,
                    )

                elapsed_samples = []
                for repetition in range(1, args.repetitions + 1):
                    elapsed, openpmd_backend = run_once()
                    elapsed_samples.append(elapsed)
                    normalized = elapsed * 1_000_000 / rays
                    row = {
                        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
                        "host": platform.node(),
                        "vtk": str(vtk_path),
                        "backend": backend,
                        "openpmd_backend": openpmd_backend,
                        "rays": rays,
                        "mode": mode,
                        "repetition": repetition,
                        "adaptive_steps": 0,
                        "rng_seed": args.rng_seed,
                        "reflection_max_iterations": args.reflection_max_iterations,
                        "reflection_tolerance": args.reflection_tolerance,
                        "elapsed_seconds": f"{elapsed:.9f}",
                        "seconds_per_million_rays": f"{normalized:.9f}",
                        "million_rays_per_second": f"{1.0 / normalized:.9f}",
                    }
                    rows.append(row)
                    print(
                        f"run {repetition}/{args.repetitions}: {backend:30s} mode={mode:10s} "
                        f"rays={rays:>9d} elapsed={elapsed:9.3f}s normalized={normalized:9.3f}s/Mray",
                        flush=True,
                    )

                median_elapsed = statistics.median(elapsed_samples)
                print(
                    f"median: {backend:30s} mode={mode:10s} rays={rays:>9d} "
                    f"elapsed={median_elapsed:9.3f}s",
                    flush=True,
                )

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
