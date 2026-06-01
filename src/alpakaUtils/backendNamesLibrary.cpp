/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <alpakaUtils/backendNames.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace hase::alpakaUtils
{
    std::vector<std::string> const& backendNames()
    {
        static std::vector<std::string> const names = hase::alpakaUtils::availableBackendNames();
        return names;
    }
} // namespace hase::alpakaUtils

extern "C" std::size_t haseAlpakaBackendCount() noexcept
{
    try
    {
        return hase::alpakaUtils::backendNames().size();
    }
    catch(...)
    {
        return 0u;
    }
}

extern "C" char const* haseAlpakaBackendName(std::size_t index) noexcept
{
    try
    {
        auto const& names = hase::alpakaUtils::backendNames();
        if(index >= names.size())
        {
            return nullptr;
        }
        return names[index].c_str();
    }
    catch(...)
    {
        return nullptr;
    }
}
