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

// STL
#include <vector>
#include <thread>

// HASEonGPU
#include <calc_phi_ase.hpp> /* calcPhiAse */
#include <types.hpp>        /* ExperimentParameters, ComputeParameters */  

unsigned calcPhiAseThreaded( const ExperimentParameters &experiment,
                         const ComputeParameters &compute,
                         Result &result){

    std::vector<std::thread> threadIds;
    
    float runtime = 0.0;

    std::vector<ComputeParameters> computes(compute.nDevices, compute);
	     
    for(unsigned device_i = 0; device_i < compute.nDevices; ++device_i){
	const unsigned samplesPerNode = compute.maxSampleRange - compute.minSampleRange+1;
	const float samplePerGpu = samplesPerNode / (float) compute.nDevices;

	computes[device_i].minSampleRange = device_i * samplePerGpu;
	computes[device_i].maxSampleRange = std::min((float)samplesPerNode, (device_i + 1) * samplePerGpu);

	threadIds.push_back(std::thread(calcPhiAse,
					std::ref(experiment),
					std::ref(compute),
					std::ref(result),
					compute.minSampleRange,
					compute.maxSampleRange,
					std::ref(runtime)));

    }

    for(unsigned i=0; i < threadIds.size(); ++i){
        threadIds[i].join();
    }
    
    return compute.nDevices;

}



