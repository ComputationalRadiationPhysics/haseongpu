#!/usr/bin/env bash

# shellcheck source=script/ci/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
set -x

versions="${OPENPMD_PROVIDER_VERSIONS:-0.17.0 0.17.1}"
work_dir="${HASE_CI_PROVIDER_WORK_DIR:-/tmp/hase-openpmd-provider-smoke}"
adios2_src="$work_dir/adios2-src"
adios2_build="$work_dir/adios2-build"
adios2_prefix="$work_dir/adios2-prefix"

rm -rf "$work_dir"
mkdir -p "$work_dir"
git clone --depth 1 --branch v2.12.1 \
    https://github.com/ornladios/ADIOS2.git "$adios2_src"
cmake -S "$adios2_src" -B "$adios2_build" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$adios2_prefix" \
    -DADIOS2_USE_MPI=OFF \
    -DADIOS2_USE_Fortran=OFF \
    -DADIOS2_USE_Python=OFF \
    -DADIOS2_USE_SST=ON \
    -DADIOS2_USE_HDF5=OFF \
    -DADIOS2_BUILD_EXAMPLES=OFF \
    -DADIOS2_BUILD_TESTING=OFF \
    -DADIOS2_INSTALL_GENERATE_CONFIG=OFF \
    -DBUILD_TESTING=OFF
cmake --build "$adios2_build" --target install --parallel "$HASE_CI_JOBS"

run_backend_smoke() {
    local openpmd_version="$1"
    local openpmd_prefix="$2"
    local backend="$3"

    ci_python "$HASE_CI_ROOT/utils/check_openpmd_compatibility.py" \
        --backend "$backend" \
        --cmake-prefix-path "$openpmd_prefix;$adios2_prefix"
    ci_python "$HASE_CI_ROOT/utils/testLaunchLaserPump.py" \
        "$backend" \
        "/tmp/hase-laser-pump-cladding-openpmd-${openpmd_version}-${backend}"
}

for openpmd_version in $versions; do
    openpmd_src="$work_dir/openpmd-src-$openpmd_version"
    openpmd_build="$work_dir/openpmd-build-$openpmd_version"
    openpmd_prefix="$work_dir/openpmd-prefix-$openpmd_version"

    git clone --depth 1 --branch "$openpmd_version" \
        https://github.com/openPMD/openPMD-api.git "$openpmd_src"
    cmake -S "$openpmd_src" -B "$openpmd_build" -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="$openpmd_prefix" \
        -DCMAKE_PREFIX_PATH="$adios2_prefix" \
        -DBUILD_TESTING=OFF \
        -DopenPMD_BUILD_EXAMPLES=OFF \
        -DopenPMD_USE_MPI=OFF \
        -DopenPMD_USE_ADIOS2=ON \
        -DopenPMD_USE_HDF5=ON \
        -DopenPMD_USE_PYTHON=OFF \
        -DADIOS2_USE_SST=ON
    cmake --build "$openpmd_build" --target install --parallel "$HASE_CI_JOBS"

    ci_pip uninstall -y HASEonGPU openpmd-api || true
    ci_pip install "openpmd-api==$openpmd_version"
    provider_cmake_args=(
        -DBUILD_TESTING=OFF
        -DHASE_TESTING=OFF
        -DDISABLE_MPI=ON
        -DHASE_OPENPMD_PROVIDER=system
        "-DHASE_RUNTIME_DIR=$work_dir/hase-runtime-$openpmd_version"
        "-DCMAKE_PREFIX_PATH=$openpmd_prefix;$adios2_prefix"
        "-DHASE_OPENPMD_RUNTIME_RPATH=$openpmd_prefix/lib;$openpmd_prefix/lib64;$adios2_prefix/lib;$adios2_prefix/lib64"
        -Dalpaka_DEP_HWLOC=ON
        -Dalpaka_HOST_MemPinningCanFail=ON
    )
    CMAKE_ARGS="${provider_cmake_args[*]}" ci_pip install "$HASE_CI_ROOT[test]"

    bash "$(dirname "${BASH_SOURCE[0]}")/smoke-installed.sh"

    mapfile -t openpmd_backends < <(ci_python "$HASE_CI_ROOT/utils/listOpenPmdBackends.py")
    backend_csv="$(ci_python "$HASE_CI_ROOT/utils/listOpenPmdBackends.py" --csv)"
    for backend in "${openpmd_backends[@]}"; do
        run_backend_smoke "$openpmd_version" "$openpmd_prefix" "$backend"
    done

    cd "$HASE_CI_ROOT"
    HASE_OPENPMD_TEST_BACKENDS="$backend_csv" ci_python -m pytest \
        tests/python/openpmd/test_openpmd_preflight.py \
        tests/python/openpmd/test_transport.py
done
