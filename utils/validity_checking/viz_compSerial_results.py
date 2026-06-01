#!/usr/bin/env python3
"""Plot serial vs Alpaka values written by compSerial_Itest."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CSV = REPO_ROOT / "tests/data/compSerial_Itest_values.csv"
DEFAULT_OUT_DIR = REPO_ROOT / "utils/validity_checking/difference_plots"


def safe_name(value: str) -> str:
    value = value.strip() or "unknown"
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def load_rows(csv_path: Path) -> dict[tuple[str, str], dict[str, list[dict[str, float | int]]]]:
    grouped: dict[tuple[str, str], dict[str, list[dict[str, float | int]]]] = defaultdict(lambda: defaultdict(list))

    with csv_path.open(newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        required = {"dataset", "backend", "vector", "index", "serial", "parallel"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{csv_path} is missing columns: {', '.join(sorted(missing))}")

        for row in reader:
            key = (row["dataset"], row["backend"])
            grouped[key][row["vector"]].append(
                {
                    "index": int(row["index"]),
                    "serial": float(row["serial"]),
                    "parallel": float(row["parallel"]),
                }
            )

    for vectors in grouped.values():
        for rows in vectors.values():
            rows.sort(key=lambda item: item["index"])

    return grouped


def plot_group(dataset: str, backend: str, vectors: dict[str, list[dict[str, float | int]]], out_dir: Path) -> Path:
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    plots = [("phiAse", axes[0]), ("dndtAse", axes[1])]
    backend_label = f"{backend or 'alpaka'} version"

    for vector_name, axis in plots:
        rows = vectors.get(vector_name, [])
        if not rows:
            axis.set_title(f"{vector_name}: no data")
            axis.grid(True, alpha=0.3)
            continue

        x = [int(row["index"]) for row in rows]
        serial = [float(row["serial"]) for row in rows]
        parallel = [float(row["parallel"]) for row in rows]

        axis.plot(x, serial, marker=".", linewidth=1.2, label="serial")
        axis.plot(x, parallel, marker=".", linewidth=1.2, label=backend_label)
        axis.set_ylabel(vector_name)
        axis.grid(True, alpha=0.3)
        axis.legend()

    axes[1].set_xlabel("index")
    figure.suptitle(f"{Path(dataset).name} - {backend or 'alpaka'}")
    figure.tight_layout()

    out_path = out_dir / f"difference_{safe_name(Path(dataset).name)}_{safe_name(backend or 'alpaka')}.png"
    figure.savefig(out_path, dpi=160)
    plt.close(figure)
    return out_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot serial and Alpaka phiAse/dndtAse values from CSV.")
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV, help=f"Input CSV; default {DEFAULT_CSV}.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR, help=f"Output directory; default {DEFAULT_OUT_DIR}.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.csv.exists():
        print(f"CSV file does not exist: {args.csv}")
        return 2

    args.out_dir.mkdir(parents=True, exist_ok=True)
    grouped = load_rows(args.csv)
    if not grouped:
        print(f"No rows found in {args.csv}")
        return 1

    for (dataset, backend), vectors in grouped.items():
        out_path = plot_group(dataset, backend, vectors, args.out_dir)
        print(f"Wrote {out_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
