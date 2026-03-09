/**
 * Copyright 2015 Erik Zenker, Carlchristian Eckert, Marius Melzer
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

#include <mesh.hpp>
#include <parser.hpp>
#include <types.hpp>
#include <simulation.hpp>
int main(int argc, char **argv){



    // Simulation data
    ExperimentParameters experiment;
    ComputeParameters    compute;
    Result               result;
    std::vector<Mesh>    meshs;

    // Parse commandline and prepate all data structures
    parse(argc, argv, experiment, compute, meshs, result);

    startSimulation<true>(experiment,compute,result,meshs);

    return 0;

}
