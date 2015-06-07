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
 * @author Marius Melzer
 * @licence GPLv3
 *
 */

#pragma once
#include <vector>

#include <mesh.hpp>
#include <types.hpp>

/**
 * @brief Calculates Phi ASE. With minRaysPerSample < maxRaysPerSample
 *        adaptive sampling can be used to improve performance.
 *
 * @param minSample_i      Smallest Index of sample point to calculate.
 * @param maxSample_i      Biggest Index of sample point to calculate.
 * @param runtime          Reference to the needed runtime.
 *
 **/
float calcPhiAse ( const ExperimentParameters& experiment,
		   const ComputeParameters& compute,
		   const Mesh& mesh,
		   Result &result,
		   const unsigned minSample_i,
		   const unsigned maxSample_i,
		   float &runtime);
