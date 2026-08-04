#!/usr/bin/env bash

# shellcheck source=script/ci/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
set -x

export DEBIAN_FRONTEND=noninteractive
apt=(apt-get)
if [[ "$(id -u)" != 0 ]]; then
    apt=(sudo apt-get)
fi
"${apt[@]}" update
"${apt[@]}" install -y python3 python3-pip python3-venv
python3 -m venv "$HASE_CI_VENV"
ci_pip install pre-commit==4.3.0
cd "$HASE_CI_ROOT"
ci_python -m pre_commit run --all-files --show-diff-on-failure
