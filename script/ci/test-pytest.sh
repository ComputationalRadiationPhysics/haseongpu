#!/usr/bin/env bash

# shellcheck source=script/ci/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
set -x

cd "$HASE_CI_ROOT"
read -r -a pytest_paths <<< "${HASE_CI_PYTEST_PATHS:-tests}"
ci_python -m pytest "${pytest_paths[@]}"
