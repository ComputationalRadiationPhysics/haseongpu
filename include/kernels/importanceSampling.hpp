/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

/**
 * @brief Forward source sampling for the volume estimator.
 *
 * No gain is estimated at sample points.  The forward estimator instead
 * chooses a tetrahedral source volume, wavelength, position, and direction
 * for each ray.
 */
#include <kernels/detail/importanceSampling.hpp>
#include <kernels/mapRaysToTets.hpp>
