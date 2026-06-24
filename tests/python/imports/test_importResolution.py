# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import json
import shutil
import subprocess
import sys
from importlib.machinery import PathFinder
from pathlib import Path

import pytest


repoRoot = Path(__file__).resolve().parents[3]


def _installedBindingsDirs():
    filtered = []
    for entry in sys.path:
        entryPath = Path(entry or ".").resolve()
        if entryPath == repoRoot:
            continue
        filtered.append(entry)

    spec = PathFinder.find_spec("HASEonGPU_Bindings", filtered)
    if spec is None:
        return []

    return [Path(location).resolve() for location in spec.submodule_search_locations or ()]


@pytest.mark.integration
def testRepoRootImportFallsBackToInstalledBindingsWithoutBuildTree(tmp_path):
    installedBindings = _installedBindingsDirs()
    if not installedBindings:
        pytest.skip("installed HASEonGPU_Bindings package is not available for this interpreter")

    checkout = tmp_path / "checkout"
    checkout.mkdir()
    shutil.copy2(repoRoot / "HASEonGPU.py", checkout / "HASEonGPU.py")
    shutil.copytree(repoRoot / "HASEonGPU_Bindings", checkout / "HASEonGPU_Bindings")
    shutil.copytree(repoRoot / "pyInclude", checkout / "pyInclude")
    runner = checkout / "smoke.py"
    runner.write_text(
        "import json\n"
        "from pathlib import Path\n"
        "import HASEonGPU\n"
        "import HASEonGPU_Bindings\n"
        "from HASEonGPU import ExperimentParameters, Result\n"
        "payload = {\n"
        "    'module': str(Path(HASEonGPU.__file__).resolve()),\n"
        "    'binding_paths': [str(Path(p).resolve()) for p in HASEonGPU_Bindings.__path__],\n"
        "    'result_module': type(Result()).__module__,\n"
        "    'min_rays': ExperimentParameters().minRaysPerSample,\n"
        "}\n"
        "print(json.dumps(payload))\n",
        encoding="utf-8",
    )

    command = [
        sys.executable,
        str(runner),
    ]
    completed = subprocess.run(
        command,
        cwd=checkout,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if completed.returncode != 0:
        pytest.fail(completed.stdout)

    payload = json.loads(completed.stdout)
    localBindings = (checkout / "HASEonGPU_Bindings").resolve()
    bindingPaths = [Path(path).resolve() for path in payload["binding_paths"]]

    assert Path(payload["module"]).resolve() == (checkout / "HASEonGPU.py").resolve()
    assert bindingPaths[0] == localBindings
    assert any(path != localBindings for path in bindingPaths)
    assert payload["result_module"] == "HASEonGPU_Bindings.HASEonGPU"
    assert payload["min_rays"] == 100000
