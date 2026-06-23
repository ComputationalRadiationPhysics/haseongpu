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


def _n_tot_by_segment(nTot, nTotGradient, topology):
    segments = int(topology.levels) - 1
    if nTotGradient is not None:
        values = np.asarray(nTotGradient, dtype=np.float64).reshape(-1)
        if values.size < segments:
            raise ValueError(f"nTotGradient must contain at least {segments} segment values")
        return values[:segments]
    if nTot is None:
        raise ValueError("calcGainFromState requires nTot or nTotGradient")
    return np.full(segments, float(nTot), dtype=np.float64)


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
    nTotGradient=None,
    roundTrip=False,
    reflectivity=1.0,
):
    """Calculate point-shaped small-signal gain directly from a ``TimeStepState``.

    The returned array has shape ``(numberOfPoints, numberOfLevels)`` and can be
    written directly as a ``vtkWedge`` point field.
    """
    topology = getattr(state, "topology", None)
    if topology is None:
        raise ValueError("calcGainFromState requires a state with topology")
    topology._require_levels()
    topology._require_thickness()

    beta_cells = np.asarray(state.betaCells, dtype=np.float64)
    expected_shape = (int(topology.numberOfPoints), int(topology.levels))
    if beta_cells.shape != expected_shape:
        if beta_cells.size != expected_shape[0] * expected_shape[1]:
            raise ValueError(f"state.betaCells must have shape {expected_shape}, got {beta_cells.shape}")
        beta_cells = beta_cells.reshape(expected_shape, order="F")

    if sigmaAbsorption is None or sigmaEmission is None:
        resolved = _resolve_cross_sections(
            crossSections=crossSections,
            spectralProperties=spectralProperties,
            laserProperties=laserProperties,
            wavelength=wavelength,
        )
        if resolved is None:
            raise ValueError("calcGainFromState requires sigmaAbsorption/sigmaEmission or crossSections")
        default_absorption, default_emission = resolved
        if sigmaAbsorption is None:
            sigmaAbsorption = default_absorption
        if sigmaEmission is None:
            sigmaEmission = default_emission

    beta_segments = 0.5 * (beta_cells[:, :-1] + beta_cells[:, 1:])
    net_absorption = float(sigmaAbsorption) - beta_segments * (float(sigmaAbsorption) + float(sigmaEmission))
    optical_depth = np.sum(
        net_absorption * _n_tot_by_segment(nTot, nTotGradient, topology)[np.newaxis, :] * float(topology.thickness),
        axis=1,
    )
    gain = np.exp(-optical_depth)
    if roundTrip:
        gain = float(reflectivity) * gain * gain
    return np.repeat(gain[:, np.newaxis], int(topology.levels), axis=1)
