/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <core/srm.hpp>

#include <charconv>
#include <cstdlib>
#include <stdexcept>
#include <string>

namespace hase::core
{
    unsigned positiveEnvironmentUnsigned(std::string_view const name, unsigned const fallback)
    {
        std::string const variable{name};
        char const* value = std::getenv(variable.c_str());
        if(value == nullptr || *value == '\0')
        {
            return fallback;
        }

        unsigned parsed = 0u;
        std::string_view const text{value};
        auto const [end, error] = std::from_chars(text.data(), text.data() + text.size(), parsed);
        if(error != std::errc{} || end != text.data() + text.size() || parsed == 0u)
        {
            throw std::runtime_error(variable + " must be a positive integer");
        }
        return parsed;
    }

    SrmControls resolveSrmControls(ExperimentParameters const& experiment)
    {
        return SrmControls{
            positiveEnvironmentUnsigned("HASE_SRM_MAX_ITERATIONS", experiment.reflectionMaxIterations),
            positiveEnvironmentUnsigned("HASE_SRM_DIVERGENCE_STREAK", 3u)};
    }

    bool srmDebugLoggingEnabled()
    {
        char const* value = std::getenv("HASE_SRM_DEBUG");
        if(value == nullptr)
        {
            return false;
        }
        std::string_view const setting{value};
        return setting == "1" || setting == "true" || setting == "TRUE" || setting == "on" || setting == "ON";
    }

    unsigned srmStatusPriority(SrmStatus const status)
    {
        switch(status)
        {
        case SrmStatus::DISABLED:
            return 0u;
        case SrmStatus::CONVERGED:
            return 1u;
        case SrmStatus::STABLE:
            return 2u;
        case SrmStatus::MAX_ITERATIONS:
            return 3u;
        case SrmStatus::DIVERGED:
            return 4u;
        }
        return 4u;
    }

} // namespace hase::core
