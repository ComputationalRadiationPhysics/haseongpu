/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 *
 * This file is part of HASENonGPU
 *
 * HASENonGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HASENonGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HASENonGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */


/**
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @author Marius Melzer
 * @licence GPLv3
 *
 */

#ifndef calc_phi_ase_H
#define calc_phi_ase_H
#include <vector>
#include <mesh.h>

/**
 * @brief Calculates Phi ASE. With minRaysPerSample < maxRaysPerSample
 *        adaptive sampling can be used to improve performance.
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
 * @param sigmaA           Vector with Absorption values
 * @param sigmaE           Vector with Emission values
 * @param mseThreshold     Threshold for adaptive and repetitive sampling.
 *                         Not reaching this threshold leads to recomputations.
 * @param useReflections   Rays can reflect on upper and lower surface of gain medium
 * @param phiAse           Reference to phiAse result (one value for every sample point).
 * @param mse              Reference to mse result (one value for every sample point).
 * @param totalRays        Reference to numberOfRays simulated per sample point.
 * @param gpu_i            Number of device that should be used.
 * @param minSample_i      Smallest Index of sample point to calculate.
 * @param maxSample_i      Biggest Index of sample point to calculate.
 * @param runtime          Reference to the needed runtime.
 *
 **/
float calcPhiAse ( const unsigned minRaysPerSample,
		   const unsigned maxRaysPerSample,
		   const unsigned maxRepetitions,
		   const Mesh& mesh,
		   const std::vector<double>& sigmaA,
		   const std::vector<double>& sigmaE,
		   const double mseThreshold,
		   const bool useReflections,
		   std::vector<float> &hPhiAse,
		   std::vector<double> &hMse,
		   std::vector<unsigned> &hTotalRays,
		   const unsigned gpu_i,
		   const unsigned minSample_i,
		   const unsigned maxSample_i,
		   float &runtime);

#endif /* calc_phi_ase_H */
