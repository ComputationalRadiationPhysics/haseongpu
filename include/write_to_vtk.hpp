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
 * @license GPLv3
 */

#pragma once
#include <vector>
#include <string>

#include <mesh.hpp>

/**
 * @brief creates a VTK file based on the mesh structure and 
 *        the dnd_ase values from a finished experiment
 *
 * @param mesh all information about triangles, points, contants. See mesh.h for details
 * @param ase the input data to write (one value for each vertex in the grid)
 * @param filename the output filename to write
 * @param minRaysPerSample with which the experiment was started (see README)
 * @param maxRaysPerSample with which the experiment was started (see README)
 * @param mseThreshold with which the experiment was started (see README)
 * @param useReflections with which the experiment was started (see README)
 * @param runtime the time needed to complete the experiment
 *
 * @return 0
 */
int writePointsToVtk(const Mesh& mesh,
	       const std::vector<double> ase,
	       const std::string filename, 
	       const unsigned minRaysPerSample,
	       const unsigned maxRaysPerSample,
	       const float mseThreshold,
	       const bool useReflections,
	       const float runtime);

/**
 * @brief creates a VTK file based on the mesh structure and 
 *        the values within prisms from a finished experiment
 *
 * @param mesh all information about triangles, points, contants. See mesh.h for details
 * @param ase the input data to write (one value for each vertex in the grid)
 * @param filename the output filename to write
 * @param minRaysPerSample with which the experiment was started (see README)
 * @param maxRaysPerSample with which the experiment was started (see README)
 * @param mseThreshold with which the experiment was started (see README)
 * @param useReflections with which the experiment was started (see README)
 * @param runtime the time needed to complete the experiment
 *
 * @return 0
 */
int writePrismToVtk(const Mesh& mesh,
	       const std::vector<double> prismData,
	       const std::string filename, 
	       const unsigned minRaysPerSample,
	       const unsigned maxRaysPerSample,
	       const float mseThreshold,
	       const bool useReflections,
	       const float runtime);

/**
 * @brief Compares a dataset with data that is extracted from a vtk-file.
 *   
 * @param compare the vector to compare with the VTK-data
 * @param filename the filename of the VTK-file to use for comparison
 * @return a vector containing the difference between "compare" and "data"
 *
 */
std::vector<double> compareVtk(std::vector<double> compare, std::string filename);



