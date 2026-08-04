#!/usr/bin/env bash

# shellcheck source=script/ci/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
set -x

cd "$HASE_CI_ROOT"
ci_python -m pytest tests
