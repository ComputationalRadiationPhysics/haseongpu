/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * HASEonGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HASEonGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HASEonGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <cstdint>
#include <limits>

namespace hase::core
{

    struct ExperimentParameters;
    class HostMesh;
    struct Result;

    class BaseVersionSerial
    {
    public:
        BaseVersionSerial(ExperimentParameters const& experiment, HostMesh& mesh, Result& result);

        void operator()(uint32_t minSampleI = 0, uint32_t maxSampleI = std::numeric_limits<uint32_t>::max());

    private:
        ExperimentParameters const& m_experiment;
        HostMesh& m_mesh;
        Result& m_result;
    };

} // namespace hase::core
