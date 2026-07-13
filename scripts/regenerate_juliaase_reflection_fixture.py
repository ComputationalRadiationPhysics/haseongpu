#!/usr/bin/env python3
"""Regenerate the committed JuliaASE reflection surface fixture.

This helper is intentionally not used by CI. It runs the repository-local
single-Tet4 driver with JuliaASE, rather than using HASE as the source of truth.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_DIR = REPO_ROOT / "tests" / "data" / "juliaASE" / "reflection_surface_reference"
REFERENCE_PATH = FIXTURE_DIR / "reference.json"
JULIA_DRIVER = REPO_ROOT / "scripts" / "juliaase_reflection_fixture.jl"


def juliaase_root() -> Path:
    configured = os.environ.get("JULIAASE_ROOT")
    root = Path(configured) if configured else REPO_ROOT.parent / "juliaASE"
    if not (root / "Project.toml").is_file():
        raise SystemExit(
            "JuliaASE checkout not found. Set JULIAASE_ROOT or place a checkout at ./juliaASE."
        )
    return root


def require_julia(root: Path) -> None:
    julia = shutil.which("julia")
    if julia is None:
        raise SystemExit("Julia executable not found on PATH; install Julia to regenerate this fixture.")
    subprocess.run(
        [julia, f"--project={root}", "-e", "using Pkg; Pkg.instantiate()"],
        check=True,
        cwd=root,
    )


def load_reference() -> dict:
    return json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))


def write_reference(reference: dict, *, dry_run: bool) -> None:
    text = json.dumps(reference, indent=2) + "\n"
    if dry_run:
        print(text)
        return
    REFERENCE_PATH.write_text(text, encoding="utf-8")


def generate_reference(root: Path, julia: str, ray_count: int) -> dict:
    environment = dict(os.environ, HASE_JULIAASE_FIXTURE_RAYS=str(ray_count))
    completed = subprocess.run(
        [julia, f"--project={root}", str(JULIA_DRIVER), str(root.resolve())],
        check=True,
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        env=environment,
    )
    try:
        return json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        raise SystemExit(f"JuliaASE fixture driver produced invalid JSON:\\n{completed.stdout}") from exc


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true", help="print regenerated JSON instead of overwriting")
    parser.add_argument(
        "--ray-count",
        type=int,
        default=1_000_000,
        help="JuliaASE histories used to establish the stored Monte Carlo reference (default: %(default)s)",
    )
    args = parser.parse_args(argv)

    if args.ray_count <= 0:
        raise SystemExit("--ray-count must be positive")

    root = juliaase_root()
    require_julia(root)
    reference = load_reference()
    generated = generate_reference(root, shutil.which("julia"), args.ray_count)
    if generated["status"] != "converged":
        raise SystemExit(f"JuliaASE fixture did not converge: {generated['status']}")

    reference["phiAse"] = generated["phiAse"]
    reference["reflectedPassWeightFractions"] = generated["reflectedPassWeightFractions"]
    reference["referenceRayCount"] = generated["rayCount"]
    reference["referenceInitialReflectedWeight"] = generated["initialReflectedWeight"]

    # JuliaASE and HASE use opposite signs for this rate convention. The stored
    # rate is deliberately the HASE/openPMD convention tested below.
    beta = np.asarray(reference["initialBetaVolume"], dtype=np.float64)
    absorption = float(max(reference["crossSections"]["crossSectionAbsorption"]))
    emission = float(max(reference["crossSections"]["crossSectionEmission"]))
    dndt = (beta * (emission + absorption) - absorption) * np.asarray(reference["phiAse"], dtype=np.float64)
    reference["dndtAse"] = dndt.tolist()
    reference["finalBetaVolume"] = (beta - float(reference["timeStep"]) * dndt).tolist()
    write_reference(reference, dry_run=args.dry_run)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
