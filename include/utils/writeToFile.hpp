/**
 * Copyright 2015 Erik Zenker, Carlchristian Eckert, Marius Melzer
 * Copyright 2026 Tim Hanel
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
 * @author Carlchristan Eckert
 * @author Erik Zenker
 * @license GPLv3
 */

#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace hase::utils
{

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
        float value,
        std::filesystem::path path,
        std::string indexName1,
        int index1,
        std::string indexName2,
        int index2);


    /**
     * @brief writes a vector to a file
     *
     * @param v the vector to write
     * @param pFilename the name of the output file
     *
     */
    void writeVectorToFile(std::vector<double> v, std::filesystem::path filename);

} // namespace hase::utils
