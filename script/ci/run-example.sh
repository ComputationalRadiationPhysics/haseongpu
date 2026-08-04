#!/usr/bin/env bash

# shellcheck source=script/ci/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
set -x

mapfile -t openpmd_backends < <(ci_python "$HASE_CI_ROOT/utils/listOpenPmdBackends.py")
backend_csv="$(ci_python "$HASE_CI_ROOT/utils/listOpenPmdBackends.py" --csv)"
ci_export HASE_OPENPMD_TEST_BACKENDS "$backend_csv"

backend="${OPENPMD_RUNTIME_BACKEND:-adios}"
if [[ ! " ${openpmd_backends[*]} " =~ " ${backend} " ]]; then
    echo "Selected runtime backend '$backend' is unavailable: ${openpmd_backends[*]}" >&2
    exit 1
fi
ci_export HASE_OPENPMD_RUNTIME_TEST_BACKEND "$backend"

ci_python "$HASE_CI_ROOT/utils/testLaunchLaserPump.py" \
    "$backend" "/tmp/hase-laser-pump-cladding-$backend"
