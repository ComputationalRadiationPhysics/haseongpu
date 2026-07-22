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
def test_installedFrontendUsesDurableOpenPmdRuntime(tmp_path):
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
        "import openpmd_api\n"
        "import pyInclude._runtime\n"
        "frontend_runtime = Path(next(iter(pyInclude._runtime.__path__))).resolve()\n"
        "runtime = pyInclude._runtime.runtime_root()\n"
        "payload = {\n"
        "  'module': str(Path(HASEonGPU.__file__).resolve()),\n"
        "  'runtime': str(runtime),\n"
        "  'frontend_vendored_calc': (frontend_runtime / 'calcPhiASE').is_file(),\n"
        "  'calc_exists': (runtime / 'calcPhiASE').is_file() or "
        "(runtime / 'python/pyInclude/_runtime/calcPhiASE').is_file(),\n"
        "  'metadata_runtime': pyInclude._runtime.runtime_config().HASE_RUNTIME_DIR,\n"
        "  'provider_dir': pyInclude._runtime.runtime_config().HASE_OPENPMD_PYTHON_PACKAGE_DIR,\n"
        "  'openpmd_module': str(Path(openpmd_api.__file__).resolve()),\n"
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
    assert not payload["frontend_vendored_calc"]
    assert Path(payload["metadata_runtime"]).resolve() == Path(payload["runtime"]).resolve()
    if payload["provider_dir"]:
        Path(payload["openpmd_module"]).resolve().relative_to(
            Path(payload["provider_dir"]).resolve()
        )
    assert not payload["legacy_bindings"]
