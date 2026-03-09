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

#include <vector>
#include <thread>

#include <mesh.hpp>
#include <calc_phi_ase.hpp>
#include <types.hpp>
#include <progressbar.hpp>

static std::vector<std::thread> threadIds;


void calcPhiAseThreaded( const ExperimentParameters &experiment,
                         const ComputeParameters &compute,
                         const Mesh& mesh,
                         Result &result,
                         const unsigned minSample_i,
                         const unsigned maxSample_i,
                         float &runtime,
                         ProgressBar &bar){

    bar.reset();
    threadIds.emplace_back(std::thread(
    [&experiment, &compute, &mesh, &result, minSample_i, maxSample_i, &runtime, &bar]()
    {
        calcPhiAse(
            experiment,
            compute,
            mesh,
            result,
            minSample_i,
            maxSample_i,
            runtime,
            bar);
    }));
}

void joinAll(){
    for(auto& t : threadIds){
        if(t.joinable()){
            t.join();
        }
    }
    threadIds.clear();
}
