/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

// Sampling a source position, spectrum, and direction is now direct forward
// volume sampling rather than a backward gain estimate from sample points.
#include <kernels/forward/volumeSampling.hpp>
