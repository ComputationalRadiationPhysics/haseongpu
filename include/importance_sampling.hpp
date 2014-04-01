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

#pragma once
#include <mesh.hpp>

/**
 * @brief Calculates preImportance which needs only to be done once.
 *
 * @param sample_i         Index of sample point for which importance sampling should be done.
 * @param reflectionSlices 1 + (2 * maxReflections) - Coded information about how many
 *                         reflections a ray should do and on which surface to start.
 * @param dMesh            All information about triangles, points, contants. 
 *                         Is located on device memory. See mesh.h for details.
 * @param sigmaA           Absorption value of the ray.
 * @param sigmaE           Emission value of the ray.
 * @param preImportance    Importance based on gain from test rays (will be returned as pointer).
 * @param blockDim         Number of threads per block.
 * @param gridDim          Number of blocks per grid.
 *
 */
void importanceSamplingPropagation(unsigned sample_i,
				   const unsigned reflectionSlices,
				   Mesh dMesh,
				   const double sigmaA,
				   const double sigmaE,
				   double *preImportance,
				   dim3 blockDim,
				   dim3 gridDim);

/**
 * @brief Calculates importance and ray distribution on prisms.
 *        Based on preImportance and triangle surfaces.
 *
 * @param reflectionSlices   1 + (2 * maxReflections) - Coded information about how many
 *                           reflections a ray should do and on which surface to start.
 * @param dMesh              All information about triangles, points, contants. 
 *                           Is located on device memory. See mesh.h for details.
 * @param raysPerSample      Number of rays that should be distributed on gain medium.
 * @param preImportance      Values calculated from importanceSamplingPropagation.
 * @param importance         Final importance values for this number of rays per sample.
 * @param raysPerPrism       Information about how many rays will start in each prism.
 * @param hSumPhi            Global memory where all threads can sum their phiAse on.
 * @param distributeRandomly In case that not all rays can be distributed by importance
 *                           sampling, the rest will be distributed randomly.
 * @param blockDim           Number of threads per block.
 * @param gridDim            Number of blocks per grid.
 *
 */
unsigned importanceSamplingDistribution(const unsigned reflectionSlices,
					Mesh dMesh,
					const unsigned raysPerSample,
					double *preImportance,
					double *importance,
					unsigned *raysPerPrism,
					float hSumPhi,
					const bool distributeRandomly,
					dim3 blockDim,
					dim3 gridDim);

