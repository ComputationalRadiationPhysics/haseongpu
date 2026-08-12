#!/usr/bin/env bash

# shellcheck source=script/ci/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
set -x

export DEBIAN_FRONTEND=noninteractive
printf 'Acquire::Retries "3";\n' > /etc/apt/apt.conf.d/80-hase-retries
apt-get update

packages=(
    ca-certificates curl git gnupg lsb-release ninja-build pkg-config
    python3 python3-dev python3-pip python3-venv
    libgl1 libglu1 libomp-dev libxcursor1 libxft2 libxinerama1 libxrender1
    software-properties-common wget
)

if [[ "${HASE_CI_HWLOC:-yes}" != no ]]; then
    packages+=(libhwloc-dev)
fi
if [[ "${HASE_CI_TBB:-off}" == on ]]; then
    packages+=(libtbb-dev)
fi
if [[ "${HASE_CI_MPI:-off}" == on ]]; then
    packages+=(openmpi-bin libopenmpi-dev)
fi
if [[ "${HASE_CI_PROVIDER_SMOKE:-0}" == 1 ]]; then
    packages+=(libhdf5-dev)
fi
if [[ "${HASE_CI_COMPILER:-gcc}" == nvcc ]]; then
    packages+=(gcc g++)
fi

apt-get install -y --no-install-recommends "${packages[@]}"
if ! command -v cmake >/dev/null 2>&1; then
    apt-get install -y --no-install-recommends cmake
fi

compiler="${HASE_CI_COMPILER:-gcc}"
version="${HASE_CI_COMPILER_VERSION:-13}"
if [[ "$compiler" == gcc ]]; then
    if command -v "gcc-$version" >/dev/null 2>&1 \
        && command -v "g++-$version" >/dev/null 2>&1; then
        :
    elif command -v gcc >/dev/null 2>&1 \
        && command -v g++ >/dev/null 2>&1 \
        && [[ "$(gcc -dumpversion)" == "$version"* ]] \
        && [[ "$(g++ -dumpversion)" == "$version"* ]]; then
        :
    else
        apt-get install -y --no-install-recommends "gcc-$version" "g++-$version"
    fi
elif [[ "$compiler" == clang && "${HASE_CI_HIP:-no}" == no ]]; then
    if ! command -v "clang-$version" >/dev/null 2>&1; then
        wget -qO /tmp/hase-llvm.sh https://apt.llvm.org/llvm.sh
        chmod +x /tmp/hase-llvm.sh
        /tmp/hase-llvm.sh "$version"
    fi
    apt-get install -y --no-install-recommends "clang-$version" "lld-$version"
fi
