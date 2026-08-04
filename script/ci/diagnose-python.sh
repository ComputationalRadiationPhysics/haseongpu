#!/usr/bin/env bash

# shellcheck source=script/ci/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
set -x

ci_python -V
ci_pip -V
ci_pip list
cd /tmp
ci_python - <<'PY'
import pathlib
import site
import sys

print("sys.executable =", sys.executable)
print("sys.path =", sys.path)
print("site-packages =", site.getsitepackages())
print("user-site =", site.getusersitepackages())
for directory in site.getsitepackages() + [site.getusersitepackages()]:
    path = pathlib.Path(directory)
    if path.exists():
        print("contents matching HASE* in", path)
        print(list(path.glob("HASE*")))
import HASEonGPU
print("HASEonGPU imported from:", HASEonGPU.__file__)
PY
