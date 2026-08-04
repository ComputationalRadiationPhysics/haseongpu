#!/usr/bin/env bash

set -euo pipefail

if [[ -z "${HASE_CI_ROOT:-}" ]]; then
    HASE_CI_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
fi

HASE_CI_BUILD_DIR="${HASE_CI_BUILD_DIR:-build/ci}"
HASE_CI_VENV="${HASE_CI_VENV:-${HASE_CI_ROOT}/.venv}"
HASE_CI_ENV_FILE="${HASE_CI_ENV_FILE:-/tmp/hase-ci-env}"
HASE_CI_JOBS="${HASE_CI_JOBS:-$(nproc)}"

if [[ -f "$HASE_CI_ENV_FILE" ]]; then
    # shellcheck disable=SC1090
    source "$HASE_CI_ENV_FILE"
fi

ci_export() {
    local name="$1"
    local value="$2"
    export "$name=$value"
    printf 'export %s=%q\n' "$name" "$value" >> "$HASE_CI_ENV_FILE"
    if [[ -n "${GITHUB_ENV:-}" ]]; then
        printf '%s=%s\n' "$name" "$value" >> "$GITHUB_ENV"
    fi
}

ci_python() {
    "$HASE_CI_VENV/bin/python" "$@"
}

ci_pip() {
    ci_python -m pip "$@"
}

ci_enabled() {
    [[ "${1:-1}" != 0 && "${1:-1}" != OFF && "${1:-1}" != off && "${1:-1}" != false ]]
}
