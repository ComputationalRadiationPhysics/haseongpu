#!/usr/bin/env bash

# This is the canonical local/container entry point for the normal CI job.
set -euo pipefail
ci_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

run_stage() {
    local toggle="$1"
    local script="$2"
    if [[ "${!toggle:-1}" != 0 ]]; then
        bash "$ci_dir/$script"
    fi
}

run_stage HASE_CI_RUN_SYSTEM_INSTALL install-system.sh
run_stage HASE_CI_RUN_SETUP setup.sh
run_stage HASE_CI_RUN_CONFIGURE configure.sh
run_stage HASE_CI_RUN_BUILD build.sh
run_stage HASE_CI_RUN_PYTHON_INSTALL install-python.sh
run_stage HASE_CI_RUN_DIAGNOSTICS diagnose-python.sh
run_stage HASE_CI_RUN_INSTALLED_SMOKE smoke-installed.sh
run_stage HASE_CI_RUN_MPI_CONFIGURE configure-mpi-tests.sh
run_stage HASE_CI_RUN_EXAMPLE run-example.sh
run_stage HASE_CI_RUN_CTEST test-ctest.sh
run_stage HASE_CI_RUN_PYTEST test-pytest.sh
