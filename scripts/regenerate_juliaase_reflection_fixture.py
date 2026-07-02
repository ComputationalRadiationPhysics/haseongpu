#!/usr/bin/env python3
"""Regenerate the committed JuliaASE reflection surface fixture.

This helper is intentionally not used by CI.  It requires a JuliaASE checkout
(via JULIAASE_ROOT or ./juliaASE) so that maintainers do not accidentally treat
HASE as the source of truth when updating the reference artifact.  The current
implementation prepares the same CI fixture schema and can be extended with the
project-local JuliaASE driver used for the reference campaign.
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


def juliaase_root() -> Path:
    configured = os.environ.get("JULIAASE_ROOT")
    root = Path(configured) if configured else REPO_ROOT / "juliaASE"
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


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true", help="print regenerated JSON instead of overwriting")
    parser.add_argument(
        "--accept-existing-reference",
        action="store_true",
        help=(
            "validate JuliaASE availability and rewrite the currently committed arrays. "
            "Use this until the local JuliaASE driver for this fixture is published."
        ),
    )
    args = parser.parse_args(argv)

    root = juliaase_root()
    require_julia(root)

    if not args.accept_existing_reference:
        raise SystemExit(
            "JuliaASE is available, but no repository-local JuliaASE driver is wired for this fixture yet. "
            "Pass --accept-existing-reference to rewrite the committed artifact unchanged, or add the "
            "JuliaASE case driver here before updating reference values."
        )

    reference = load_reference()
    # Recompute the stored one-step beta update from the reference dndt to keep
    # deterministic formatting and catch accidental manual edits.
    beta = np.asarray(reference["initialBetaVolume"], dtype=np.float64)
    dndt = np.asarray(reference["dndtAse"], dtype=np.float64)
    reference["finalBetaVolume"] = (beta - float(reference["timeStep"]) * dndt).tolist()
    write_reference(reference, dry_run=args.dry_run)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
