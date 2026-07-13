/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <core/types.hpp>

#include <string_view>

namespace hase::core
{
    struct SrmControls
    {
        unsigned maxIterations;
        unsigned divergenceStreak;
    };

    [[nodiscard]] unsigned positiveEnvironmentUnsigned(std::string_view name, unsigned fallback);

    [[nodiscard]] SrmControls resolveSrmControls(ExperimentParameters const& experiment);

    [[nodiscard]] bool srmDebugLoggingEnabled();

    [[nodiscard]] unsigned srmStatusPriority(SrmStatus status);

} // namespace hase::core
