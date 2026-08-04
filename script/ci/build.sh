#!/usr/bin/env bash

# shellcheck source=script/ci/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
set -x

cmake --build "$HASE_CI_ROOT/$HASE_CI_BUILD_DIR" --parallel "$HASE_CI_JOBS"
