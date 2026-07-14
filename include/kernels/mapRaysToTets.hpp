/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

/**
 * @brief Directly maps each forward ray to a sampled tetrahedron.
 *
 * The former prism-distribution algorithm is not used by the forward volume
 * estimator.  Its replacement samples tetrahedra by active volume or by
 * beta-weighted volume before launching a ray.
 */
#include <kernels/detail/mapRaysToTets.hpp>
