/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

// Forward rays are sampled directly in tetrahedral volumes.  Keep the
// historical mapping layer at a tetrahedron-specific path while the actual
// device helpers remain grouped with the forward implementation.
#include <kernels/forward/volumeSampling.hpp>
