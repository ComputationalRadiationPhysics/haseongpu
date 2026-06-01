/**
 * Copyright 2015 Erik Zenker, Carlchristian Eckert, Marius Melzer
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


/**
 * @mainpage HASEonGPU - High performance Amplified Spontaneous EmissioN on GPU
 *
 * Project with HZDR for porting their ASE-code to a GPU cluster.
 *
 * @author Erik Zenker, Carlchristian Eckert, Marius Melzer
 */

// STL

#include <alpaka/alpaka.hpp>

#include <core/mesh.hpp>
#include <core/simulation.hpp>
#include <core/types.hpp>
#include <parse/parser.hpp>

#include <exception>
#include <iostream>

int main(int argc, char** argv)
{
    try
    {
        // Simulation data
        hase::core::ExperimentParameters experiment;
        hase::core::ComputeParameters compute;
        hase::core::Result result;
        std::optional<hase::core::HostMesh> mesh;
        // Parse commandline and prepate all data structures
        hase::parse::parse(argc, argv, experiment, compute, mesh, result);
        return hase::core::startSimulation<true>(experiment, compute, result, *mesh);
    }
    catch(std::exception const& e)
    {
        std::cerr << "[ERROR] " << e.what() << std::endl;
        return 1;
    }
    catch(...)
    {
        std::cerr << "[ERROR] Unknown exception" << std::endl;
        return 1;
    }


    return 0;
}
