/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

/**
 * @brief Tetrahedral forward-ray traversal helpers.
 *
 * Rays advance face by face through the explicit Tet4 mesh and accumulate
 * their contribution in the forward estimator.  The device functions stay
 * header-only so every Alpaka backend can compile them for its accelerator.
 */
#include <kernels/forward/rayWalk.hpp>
