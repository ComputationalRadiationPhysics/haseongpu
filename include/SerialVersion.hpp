/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
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

#include <mesh.hpp>
#include <types.hpp>

#include <vector>
#pragma once

// only forward definitions
namespace
{
    struct BaseVersionSerialContext;
    struct PointM;
} // namespace
struct ExperimentParameters;
class HostMesh;
struct Result;

class BaseVersionSerial
{
public:
    BaseVersionSerial(ExperimentParameters const& experiment, HostMesh& mesh, Result& result);

    void operator()(uint32_t minSampleI = 0, uint32_t maxSampleI = std::numeric_limits<uint32_t>::max());

private:
    void mainLoop(
        BaseVersionSerialContext& ctx,
        std::vector<PointM> const& points,
        std::vector<double>& dndtAse,
        std::vector<double>& betaCells,
        float hostCrystalFluorescence,
        uint32_t minSampleI,
        uint32_t maxSampleI);
    ExperimentParameters const& m_experiment;
    HostMesh& m_mesh;
    Result& m_result;
};
