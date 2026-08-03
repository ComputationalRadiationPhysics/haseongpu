# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np
import pytest

from HASEonGPU import TimeStepState, calcGainFromState


def _state(topology, beta):
    return TimeStepState(
        step=1,
        time=1e-9,
        betaVolume=np.asarray(beta, dtype=np.float64),
        phiAse=None,
        standardError=None,
        relativeStandardError=None,
        totalRays=None,
        dndtAse=np.zeros((topology.numberOfTriangles, topology.levels - 1)),
        dndtPump=np.zeros((topology.numberOfTriangles, topology.levels - 1)),
        aseResult=None,
        topology=topology,
    )


def test_calcGainFromStateReturnsVtkCellField(smallTopology, crossSections):
    beta = np.full((smallTopology.numberOfTriangles, smallTopology.levels - 1), 0.25)
    state = _state(smallTopology, beta)

    gain = calcGainFromState(state, crossSections, nTot=2.0)

    sigma_abs = crossSections.crossSectionAbsorption[0]
    sigma_ems = crossSections.crossSectionEmission[0]
    expected = (0.25 * (sigma_abs + sigma_ems) - sigma_abs) * 2.0
    assert gain.shape == (smallTopology.numberOfPrisms,)
    assert np.allclose(gain, expected)


def test_calcGainFromStateRequiresNTot(smallTopology, crossSections):
    beta = np.ones((smallTopology.numberOfTriangles, smallTopology.levels - 1), dtype=np.float64)
    state = _state(smallTopology, beta)

    with pytest.raises(ValueError, match="requires nTot"):
        calcGainFromState(state, crossSections)
