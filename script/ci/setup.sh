#!/usr/bin/env bash

# shellcheck source=script/ci/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
set -x

: > "$HASE_CI_ENV_FILE"
python3 -m venv "$HASE_CI_VENV"
ci_pip install --upgrade pip

compiler="${HASE_CI_COMPILER:-gcc}"
version="${HASE_CI_COMPILER_VERSION:-13}"
cc=
cxx=
cudacxx=

if [[ "$compiler" == gcc ]]; then
    cc="$(command -v "gcc-$version" 2>/dev/null || command -v gcc)"
    cxx="$(command -v "g++-$version" 2>/dev/null || command -v g++)"
elif [[ "$compiler" == clang ]]; then
    if [[ "${HASE_CI_HIP:-no}" != no && -x /opt/rocm/llvm/bin/clang ]]; then
        cc=/opt/rocm/llvm/bin/clang
        cxx=/opt/rocm/llvm/bin/clang++
    else
        cc="$(command -v "clang-$version" 2>/dev/null || command -v clang)"
        cxx="$(command -v "clang++-$version" 2>/dev/null || command -v clang++)"
    fi
    if [[ "${HASE_CI_CUDA:-no}" != no ]]; then
        cudacxx="$cxx"
    fi
elif [[ "$compiler" == nvcc ]]; then
    cc="$(command -v gcc)"
    cxx="$(command -v g++)"
    cudacxx="$(command -v nvcc)"
else
    echo "Unsupported HASE_CI_COMPILER=$compiler" >&2
    exit 1
fi

if [[ "$compiler" != nvcc && "$($cc -dumpversion)" != "$version"* ]]; then
    echo "Selected compiler $cc does not provide requested version $version" >&2
    exit 1
fi

ci_export CC "$cc"
ci_export CXX "$cxx"
if [[ -n "$cudacxx" ]]; then
    ci_export CUDACXX "$cudacxx"
fi

if [[ "${HASE_CI_HIP:-no}" != no ]]; then
    ci_export ROCM_PATH "${ROCM_PATH:-/opt/rocm}"
    if [[ -x /opt/rocm/llvm/bin/clang++ ]]; then
        ci_export HIPCXX /opt/rocm/llvm/bin/clang++
    else
        ci_export HIPCXX /opt/rocm/bin/hipcc
    fi
fi

runtime_paths=()
for cuda_dir in /usr/local/cuda /usr/local/cuda-*; do
    if [[ -d "$cuda_dir/lib64" ]]; then
        runtime_paths+=("$cuda_dir/lib64")
    fi
done
if [[ -d /opt/rocm/lib ]]; then
    runtime_paths+=(/opt/rocm/lib)
fi
if (( ${#runtime_paths[@]} > 0 )); then
    joined="$(IFS=:; echo "${runtime_paths[*]}")"
    ci_export LD_LIBRARY_PATH "${joined}:${LD_LIBRARY_PATH:-}"
fi

if [[ "${HASE_CI_MPI:-off}" == on ]]; then
    ci_export OMPI_ALLOW_RUN_AS_ROOT 1
    ci_export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM 1
    ci_export HASE_MPIEXEC_EXTRA_ARGS --oversubscribe
fi

if [[ -n "${GITHUB_PATH:-}" ]]; then
    printf '%s\n' "$HASE_CI_VENV/bin" >> "$GITHUB_PATH"
fi
