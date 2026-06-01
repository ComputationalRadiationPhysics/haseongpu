# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later


import numpy as np
import pytest

from HASEonGPU import LaserProperties


def testLaserPropertiesDescribeAndExportCalcPhiAseDict():
    laser = LaserProperties.spectral(
        wavelengthsAbsorption=[900.0, 910.0],
        crossSectionAbsorption=[0.01, 0.02],
        wavelengthsEmission=[1000.0, 1010.0],
        crossSectionEmission=[0.03, 0.04],
        resolution=2,
    )

    sigmaAbs = laser.get("sigmaA")
    spectralResolution = laser.get("spectralResolution")

    assert sigmaAbs.name == "s_abs"
    assert sigmaAbs.expectedShape == ("nAbsorptionSamples",)
    assert "Absorption cross-section" in sigmaAbs.description
    assert spectralResolution.value == 2
    assert laser.maxSigmaA == 0.02
    assert laser.maxSigmaE == 0.04

    exported = laser.toDict()

    assert set(exported) == {"l_abs", "l_ems", "s_abs", "s_ems", "l_res"}
    assert exported["l_res"] == 2
    assert np.allclose(exported["s_ems"], [0.03, 0.04])


def testLaserPropertiesValidateMatchingSpectrumLengths():
    laser = LaserProperties()
    laser.get("l_abs").value = [900.0, 910.0]

    with pytest.raises(ValueError, match="l_abs and s_abs"):
        laser.get("s_abs").value = [0.01]


def testLaserPropertiesMonochromaticConstructor():
    laser = LaserProperties.monochromatic(absorption=0.01, emission=0.02)

    assert laser.toDict()["l_res"] == 1
    assert np.allclose(laser.toDict()["s_abs"], [0.01])
    assert np.allclose(laser.toDict()["s_ems"], [0.02])
