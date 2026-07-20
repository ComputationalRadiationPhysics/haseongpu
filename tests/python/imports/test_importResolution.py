# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import json
import os
import subprocess
import sys
from importlib.metadata import PackageNotFoundError, distribution
from pathlib import Path

import pytest


repoRoot = Path(__file__).resolve().parents[3]


@pytest.mark.integration
def test_installedFrontendHasPrivateOpenPmdRuntime(tmp_path):
    try:
        distribution("HASEonGPU")
    except PackageNotFoundError:
        pytest.skip("HASEonGPU wheel is not installed for this interpreter")

    runner = tmp_path / "smoke.py"
    runner.write_text(
        "import importlib.util\n"
        "import json\n"
        "from pathlib import Path\n"
        "import HASEonGPU\n"
        "import pyInclude._runtime\n"
        "runtime = Path(next(iter(pyInclude._runtime.__path__))).resolve()\n"
        "payload = {\n"
        "  'module': str(Path(HASEonGPU.__file__).resolve()),\n"
        "  'runtime': str(runtime),\n"
        "  'calc_exists': (runtime / 'calcPhiASE').is_file(),\n"
        "  'legacy_bindings': importlib.util.find_spec('HASEonGPU_Bindings') is not None,\n"
        "}\n"
        "print(json.dumps(payload))\n",
        encoding="utf-8",
    )
    env = os.environ.copy()
    env["PYTHONNOUSERSITE"] = "1"
    env.pop("PYTHONPATH", None)
    completed = subprocess.run(
        [sys.executable, str(runner)],
        cwd=tmp_path,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if completed.returncode != 0:
        pytest.fail(completed.stdout)

    payload = json.loads(completed.stdout)
    assert "site-packages" in payload["module"] or "dist-packages" in payload["module"]
    assert payload["calc_exists"]
    assert not payload["legacy_bindings"]
