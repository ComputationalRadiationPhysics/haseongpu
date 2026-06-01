# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

from pathlib import Path

_build_bindings = Path(__file__).resolve().parents[1] / "build" / "python" / "HASEonGPU_Bindings"
if _build_bindings.is_dir():
    __path__.insert(0, str(_build_bindings))

from .HASEonGPU import *
