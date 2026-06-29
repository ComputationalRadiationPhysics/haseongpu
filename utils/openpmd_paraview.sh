#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage:
  utils/openpmd_paraview.sh <openPMD .pmd/.h5 file> [paraview args...]

Environment overrides:
  HASE_REPO_ROOT
  HASE_PARAVIEW_OPENPMD_PYTHONPATH
  HASE_PARAVIEW_HDF5_LIB

The helper disables the Python user site to avoid mixing user-installed VTK
with ParaView's bundled VTK, then injects the HASE-built HDF5-enabled
openpmd_api and matching HDF5 library.
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -lt 1 ]]; then
    usage
    exit 0
fi

target_file="$1"
shift

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="${HASE_REPO_ROOT:-$(cd -- "${script_dir}/.." && pwd)}"
openpmd_pythonpath="${HASE_PARAVIEW_OPENPMD_PYTHONPATH:-${repo_root}/build/codex-openpmd-hdf5-superbuild/_deps/openpmd-build/lib/python3.12/site-packages}"
hdf5_lib="${HASE_PARAVIEW_HDF5_LIB:-${repo_root}/build/codex-openpmd-hdf5-superbuild/_deps/hdf5-build/src/libhdf5.so.310}"

if [[ ! -e "${target_file}" ]]; then
    echo "openPMD file does not exist: ${target_file}" >&2
    exit 2
fi
if [[ ! -d "${openpmd_pythonpath}" ]]; then
    echo "openpmd_api Python path does not exist: ${openpmd_pythonpath}" >&2
    exit 2
fi
if [[ ! -e "${hdf5_lib}" ]]; then
    echo "HDF5 library does not exist: ${hdf5_lib}" >&2
    exit 2
fi

export PYTHONNOUSERSITE=1
export PYTHONPATH="${openpmd_pythonpath}${PYTHONPATH:+:${PYTHONPATH}}"
export LD_PRELOAD="${hdf5_lib}${LD_PRELOAD:+:${LD_PRELOAD}}"

exec paraview "${target_file}" "$@"
