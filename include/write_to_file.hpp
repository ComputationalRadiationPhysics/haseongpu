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
 * @author Carlchristan Eckert
 * @author Erik Zenker 
 * @license GPLv3
 */

#pragma once
#include <string>
#include <vector>

/** 
 * @brief writes a value to a file, where the filename is appended with 2 longish indices
 *
 * @param value the value to write to a file
 * @param path the path where to create the file
 * @param indexName1 identifier of the first index
 * @param index1 the value of the first index
 * @param indexName2 identifier of the second index
 * @param index2 the value of the second index
 */
int writeValueToFile(
    const float value, 
    const std::string path, 
    const std::string indexName1, 
    const int index1, 
    const std::string indexName2, 
    const int index2
    );


/**
 * @brief writes a vector to a file
 *
 * @param v the vector to write
 * @param pFilename the name of the output file
 *
 */
void writeVectorToFile(std::vector<double> v, std::string pFilename);
