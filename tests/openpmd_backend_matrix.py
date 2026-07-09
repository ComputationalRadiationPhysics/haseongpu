# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import os

BACKEND_PRIORITY = ("adios-sst", "adios", "hdf5")


def _clean(values):
    selected = {str(value).strip().lower() for value in values if str(value).strip()}
    return [backend for backend in BACKEND_PRIORITY if backend in selected]


def _env_backends():
    return _clean(os.environ.get("HASE_OPENPMD_TEST_BACKENDS", "").split(","))


def _front_end_backends():
    from HASEonGPU import OpenPmdBackends

    return OpenPmdBackends.all()


def openpmd_test_backends():
    backends = _env_backends() or _clean(_front_end_backends())
    if not backends:
        raise RuntimeError(
            "The HASEonGPU frontend did not report any available openPMD backends. "
            "Build HaseOpenPmdBackendProbe or set HASE_OPENPMD_TEST_BACKENDS for an explicit manual matrix."
        )
    return backends


def openpmd_runtime_backend():
    return openpmd_test_backends()[0]


def openpmd_runtime_executable():
    from pyInclude.openpmd.transport import findCalcPhiAse

    return findCalcPhiAse()
