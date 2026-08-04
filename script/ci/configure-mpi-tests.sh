#!/usr/bin/env bash

# shellcheck source=script/ci/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
set -x

if [[ "${HASE_CI_MPI:-off}" != on || -z "${HASE_CI_MPI_RANKS:-}" ]]; then
    exit 0
fi

config_dir="$HASE_CI_ROOT/$HASE_CI_BUILD_DIR/ci-config"
mkdir -p "$config_dir"
ci_python - "$HASE_CI_ROOT" "$config_dir" "$HASE_CI_MPI_RANKS" <<'PY'
import sys
from pathlib import Path

import yaml

root = Path(sys.argv[1])
destination = Path(sys.argv[2])
ranks = int(sys.argv[3])
for relative in (
    Path("config/hase-phiase-mpi.yaml"),
    Path("tests/data/cfg/phiAseTestConfig-mpi.yaml"),
):
    config = yaml.safe_load((root / relative).read_text(encoding="utf-8"))
    simulation = config.get("simulation")
    if not isinstance(simulation, dict):
        raise RuntimeError(f"{relative} does not contain a simulation mapping")
    phi_ase = simulation.get("phi_ase")
    if not isinstance(phi_ase, dict):
        raise RuntimeError(f"{relative} does not contain simulation.phi_ase")
    phi_ase["parallel_mode"] = "mpi"
    phi_ase["n_per_node"] = ranks
    output = destination / relative.name
    output.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    print(f"{output}: parallel_mode=mpi, n_per_node={ranks}")
PY

ci_export HASE_PHIASE_CONFIG "$config_dir/hase-phiase-mpi.yaml"
ci_export HASE_TEST_PHIASE_CONFIG "$config_dir/phiAseTestConfig-mpi.yaml"
