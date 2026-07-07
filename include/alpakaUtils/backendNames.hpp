/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <alpaka/alpaka.hpp>

#include <string>
#include <vector>

namespace hase::alpakaUtils
{
    inline std::string getNameForBackend(auto const& backend, auto const& device)
    {
        std::string backendName;
        backendName += alpaka::onHost::getName(alpaka::getApi(device)) + "_";
        backendName += alpaka::onHost::getName(alpaka::getDeviceKind(device)) + "_";
        backendName += alpaka::onHost::getName(backend[alpaka::object::exec]);
        return backendName;
    }

    inline std::vector<std::string> availableBackendNames()
    {
        auto backends = alpaka::onHost::allBackends(alpaka::onHost::enabledApis, alpaka::exec::enabledExecutors);
        std::vector<std::string> names;
        alpaka::onHost::executeForEachIfHasDevice(
            [&](auto const& backend) -> int
            {
                auto devSelector = alpaka::onHost::makeDeviceSelector(backend[alpaka::object::deviceSpec]);
                auto sampleDevice = devSelector.makeDevice(0);
                names.emplace_back(getNameForBackend(backend, sampleDevice));
                return 0;
            },
            backends);
        return names;
    }
} // namespace hase::alpakaUtils
