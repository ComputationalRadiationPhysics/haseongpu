# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Gain calculation helpers for simulation state snapshots."""

from __future__ import annotations

import numpy as np

from .laser import CrossSectionData, LaserProperties


def _cross_sections_from_laser_properties(laser_properties):
    values = laser_properties.toDict()
    return CrossSectionData(
        wavelengthsAbsorption=values["l_abs"],
        crossSectionAbsorption=values["s_abs"],
        wavelengthsEmission=values["l_ems"],
        crossSectionEmission=values["s_ems"],
        resolution=values["l_res"],
    )


def _resolve_cross_sections(crossSections=None, spectralProperties=None, laserProperties=None, wavelength=None):
    spectra = crossSections if crossSections is not None else spectralProperties
    if spectra is None and laserProperties is not None:
        spectra = _cross_sections_from_laser_properties(laserProperties)
    if isinstance(spectra, LaserProperties):
        spectra = _cross_sections_from_laser_properties(spectra)
    if spectra is None:
        return None
    if not isinstance(spectra, CrossSectionData):
        raise TypeError("crossSections must be CrossSectionData, SpectralDecomposition, or LaserProperties")

    if wavelength is None:
        peak_index = int(np.argmax(spectra.crossSectionEmission))
        wavelength = float(spectra.wavelengthsEmission[peak_index])
    return spectra.absorptionAt(wavelength), spectra.emissionAt(wavelength)


def calcGainFromState(
    state,
    crossSections=None,
    nTot=None,
    *,
    spectralProperties=None,
    laserProperties=None,
    wavelength=None,
    sigmaAbsorption=None,
    sigmaEmission=None,
):
    """Calculate point-shaped local small-signal gain from a ``TimeStepState``.

    The returned array has shape ``(numberOfPoints, numberOfLevels)`` and can be
    written directly as a ``vtkWedge`` point field.
    """
    topology = getattr(state, "topology", None)
    if topology is None:
        raise ValueError("calcGainFromState requires a state with topology")
    topology._require_levels()

    beta_cells = np.asarray(state.betaCells, dtype=np.float64)
    expected_shape = (int(topology.numberOfPoints), int(topology.levels))
    if beta_cells.shape != expected_shape:
        if beta_cells.size != expected_shape[0] * expected_shape[1]:
            raise ValueError(f"state.betaCells must have shape {expected_shape}, got {beta_cells.shape}")
        beta_cells = beta_cells.reshape(expected_shape, order="F")

    if nTot is None:
        raise ValueError("calcGainFromState requires nTot for local gain")

    sigma_absorption = sigmaAbsorption
    sigma_emission = sigmaEmission
    if sigma_absorption is None or sigma_emission is None:
        resolved = _resolve_cross_sections(
            crossSections=crossSections,
            spectralProperties=spectralProperties,
            laserProperties=laserProperties,
            wavelength=wavelength,
        )
        if resolved is None:
            raise ValueError("calcGainFromState requires sigmaAbsorption/sigmaEmission or crossSections")
        default_absorption, default_emission = resolved
        if sigma_absorption is None:
            sigma_absorption = default_absorption
        if sigma_emission is None:
            sigma_emission = default_emission

    local_gain = (
        beta_cells * (float(sigma_absorption) + float(sigma_emission)) - float(sigma_absorption)
    ) * float(nTot)
    return local_gain
