#!/usr/bin/env bash

# shellcheck source=script/ci/cmake-args.sh
source "$(dirname "${BASH_SOURCE[0]}")/cmake-args.sh"
set -x

make_hase_cmake_args
cmake -S "$HASE_CI_ROOT" -B "$HASE_CI_ROOT/$HASE_CI_BUILD_DIR" -G Ninja \
    "${HASE_CMAKE_ARGS[@]}" \
    -DHASE_OPENPMD_PROVIDER=bundled
