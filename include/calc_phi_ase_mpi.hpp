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


/**
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @licence GPLv3
 *
 */

#pragma once
#include <vector>

#include <mesh.hpp>

/**
 * @brief A wrapper for calcPhiAse, that distributes sample points
 *        to the available MPI nodes.The Nodes will split 
 *        up in one head node and the others as compute nodes. 
 *        The head node distributes the available sample
 *        points by demand.
 *
 * @param minRaysPerSample Lower bound for raysPerSample
 *                         in case of adaptive sampling.
 * @param maxRaysPerSample Uppper boud for raysPerSample
 *                         in case of adaptive sampling.
 * @param maxRepetitions   Number of Repetitions will
 *                         be done, when not reaching mse threshold
 * @param dMesh            All information about triangles, points, contants. 
 *                         Is located in device memory. See mesh.h for details.
 * @param hMesh            Same as dMesh, but locatet in host memory.
 * @param hsigmaA           Vector with Absorption values
 * @param hsigmaE           Vector with Emission values
 * @param mseThreshold     Threshold for adaptive and repetitive sampling.
 *                         Not reaching this threshold leads to recomputations.
 * @param useReflections   Rays can reflect on upper and lower surface of gain medium
 * @param phiAse           Reference to phiAse result (one value for every sample point).
 * @param mse              Reference to mse result (one value for every sample point).
 * @param totalRays        Reference to numberOfRays simulated per sample point.
 * @param gpu_i            Number of device that should be used.
 *
 * @return number of used compute nodes
 */
float calcPhiAseMPI ( const unsigned minRaysPerSample,
		      const unsigned maxRaysPerSample,
		      const unsigned maxRepetitions,
		      const Mesh& mesh,
		      const std::vector<double>& hSigmaA,
		      const std::vector<double>& hSigmaE,
		      const double mseThreshold,
		      const bool useReflections,
		      std::vector<float> &hPhiAse,
		      std::vector<double> &hMse,
		      std::vector<unsigned> &hTotalRays,
		      const unsigned gpu_i);

