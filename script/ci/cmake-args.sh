#!/usr/bin/env bash

# shellcheck source=script/ci/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

make_hase_cmake_args() {
    HASE_CMAKE_ARGS=(
        -DBUILD_TESTING=ON
        -DHASE_TESTING=ON
    )

    if [[ "${HASE_CI_MPI:-off}" == on ]]; then
        HASE_CMAKE_ARGS+=(-DDISABLE_MPI=OFF)
    else
        HASE_CMAKE_ARGS+=(-DDISABLE_MPI=ON)
    fi

    case "${OPENPMD_RUNTIME_BACKEND:-adios}" in
        hdf5)
            HASE_CMAKE_ARGS+=(
                -DHASE_OPENPMD_USE_ADIOS2=OFF
                -DHASE_OPENPMD_USE_SST=OFF
                -DHASE_OPENPMD_USE_HDF5=ON
                -DHASE_OPENPMD_FETCH_HDF5=ON
            )
            ;;
        adios-sst)
            HASE_CMAKE_ARGS+=(
                -DHASE_OPENPMD_USE_ADIOS2=ON
                -DHASE_OPENPMD_USE_SST=ON
                -DHASE_OPENPMD_USE_HDF5=OFF
            )
            ;;
        adios)
            HASE_CMAKE_ARGS+=(
                -DHASE_OPENPMD_USE_ADIOS2=ON
                -DHASE_OPENPMD_USE_SST=OFF
                -DHASE_OPENPMD_USE_HDF5=OFF
            )
            ;;
        *)
            echo "Unsupported OPENPMD_RUNTIME_BACKEND=${OPENPMD_RUNTIME_BACKEND}" >&2
            return 1
            ;;
    esac

    if [[ "${HASE_CI_HWLOC:-yes}" == no ]]; then
        HASE_CMAKE_ARGS+=(-Dalpaka_DEP_HWLOC=OFF)
    else
        HASE_CMAKE_ARGS+=(
            -Dalpaka_DEP_HWLOC=ON
            -Dalpaka_HOST_MemPinningCanFail=ON
        )
    fi

    if [[ "${HASE_CI_SIMD:-DEFAULT}" == EMULATION ]]; then
        HASE_CMAKE_ARGS+=(-Dalpaka_SIMD=EMULATION)
    fi

    if [[ "${HASE_CI_HIP:-no}" != no ]]; then
        HASE_CMAKE_ARGS+=(
            "-DCMAKE_PREFIX_PATH=${ROCM_PATH:-/opt/rocm}"
            "-DCMAKE_HIP_COMPILER=${HIPCXX:-/opt/rocm/llvm/bin/clang++}"
            -DCMAKE_HIP_ARCHITECTURES=gfx900
        )
    fi
}
