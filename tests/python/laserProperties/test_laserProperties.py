# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later


import numpy as np
import pytest

from HASEonGPU import CrossSectionData, LaserProperties


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


def testCrossSectionDataExposesSpectralFields():
    cross_sections = CrossSectionData(
        wavelengthsAbsorption=[900.0, 910.0],
        crossSectionAbsorption=[0.01, 0.02],
        wavelengthsEmission=[1000.0, 1010.0],
        crossSectionEmission=[0.03, 0.04],
        resolution=2,
    )

    fields = {field.name: field for field in cross_sections.getFields()}

    assert set(fields) == {"lambdaAbsorption", "lambdaEmission", "sigmaAbsorption", "sigmaEmission"}
    assert fields["lambdaAbsorption"].meta()["recordName"] == "lambda_absorption"
    assert fields["lambdaAbsorption"].meta()["unit"] == "m"
    assert fields["sigmaAbsorption"].meta()["unit"] == "cm^2"
    assert fields["sigmaAbsorption"].value().tolist() == [0.01, 0.02]

    cross_sections.getField("sigmaA").value([0.05, 0.06])
    np.testing.assert_array_equal(cross_sections.crossSectionAbsorption, [0.05, 0.06])

    with pytest.raises(ValueError, match="wavelengthsAbsorption and crossSectionAbsorption"):
        cross_sections.getField("sigmaAbsorption").value([0.07])


def testCrossSectionDataOpenPmdFieldsAreFieldObjects():
    cross_sections = CrossSectionData.monochromatic(
        wavelength=940e-9,
        crossSectionAbsorption=0.01,
        crossSectionEmission=0.02,
    )

    fields = list(cross_sections.openPmdFields(lambda values: type("Context", (), {"spectral": len(values)})()))

    assert [field.name for field in fields] == [
        "lambdaAbsorption",
        "lambdaEmission",
        "sigmaAbsorption",
        "sigmaEmission",
    ]
    assert fields[0].spec.recordName == "lambda_absorption"
    assert fields[2].spec.recordName == "sigma_absorption"
    assert fields[0].context.spectral == 1
