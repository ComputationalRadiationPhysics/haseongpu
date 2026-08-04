#!/usr/bin/env bash

# shellcheck source=script/ci/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
set -x

cp "$HASE_CI_ROOT/tests/python/imports/hase_runtime_smoketest.py" /tmp/hase_runtime_smoketest.py
cp "$HASE_CI_ROOT/tests/python/imports/hase_frontend_smoketest.py" /tmp/hase_frontend_smoketest.py
cd /tmp
ci_python hase_runtime_smoketest.py
ci_python hase_frontend_smoketest.py
