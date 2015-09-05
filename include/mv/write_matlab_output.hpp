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


#pragma once

#include <vector>

#include <boost/filesystem/path.hpp> /* boost::filesystem::path */

/**
 * @brief creates textfiles containing all the results as a 2D matrix
 *        (needs to be reshaped in matlab, see calcPhiASE.m)
 *
 * @param experimentPath the path where to create the output files
 * @param ase the phi_ASE values to print
 * @param N_rays the number of rays used per sample point
 * @param mse_values the error for each sample point
 * @param numberOfSamples the amount of vertices that were sampled in the mesh (over all levels)
 * @param numberOfLevels the amount of levels in which the gain medium was split 
 *
 * @author Carlchristian Eckert
 * @author Erik Zenker
 * @author Marius Melzer
 *
 * @license GPLv3
 */
void writeMatlabOutput(
    const boost::filesystem::path experimentPath,
    const std::vector<float> ase,
    const std::vector<unsigned> N_rays, 
    const std::vector<double> mse_values,
    const unsigned numberOfSamples,
    const unsigned numberOfLevels
    );

