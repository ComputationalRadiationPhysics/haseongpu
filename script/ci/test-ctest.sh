#!/usr/bin/env bash

# shellcheck source=script/ci/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
set -x

ctest --test-dir "$HASE_CI_ROOT/$HASE_CI_BUILD_DIR" --output-on-failure
