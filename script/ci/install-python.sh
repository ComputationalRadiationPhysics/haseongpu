#!/usr/bin/env bash

# shellcheck source=script/ci/cmake-args.sh
source "$(dirname "${BASH_SOURCE[0]}")/cmake-args.sh"
set -x

make_hase_cmake_args
HASE_CMAKE_ARGS+=(
    -DHASE_BUILD_OPENPMD_FROM_SOURCE=ON
    "-DHASE_RUNTIME_DIR=$HASE_CI_ROOT/$HASE_CI_BUILD_DIR"
)

CMAKE_ARGS="${HASE_CMAKE_ARGS[*]}" ci_pip install -v "$HASE_CI_ROOT[test]"
