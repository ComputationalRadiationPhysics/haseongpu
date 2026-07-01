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


#include <utils/writeMatlabOutput.hpp>

#include <fstream>
#include <iomanip> /* std::setprecision(), std::fixed() */
#include <vector>

namespace hase::utils
{

    /**
     * @brief write input data in a 3D-matrix format compatible with reshape-style post-processing
     *
     * @param data the data to write into the matrix (as 1D vector)
     *        layout of the vector: (2D matrix for z==0 in row-major order), (... for z==1 in row-major order) , ...
     * @param file the destination file for the matrix
     * @param rowCount the number of rows for the output matrix (x-dimension)
     * @param columnCount the number of columns for the output matrix (y-dimension)
     * @param pageCount the number of pages for the output matrix (z-dimension)
     *
     */
    template<typename T>
    void write3dMatrix(
        std::vector<T> const& data,
        std::ofstream& file,
        unsigned const rowCount,
        unsigned const columnCount,
        unsigned const pageCount)
    {
        unsigned elementsPerPage = rowCount * columnCount;
        // write first line, containing geometry information parsable by matlab
        file << rowCount << " " << columnCount << " " << pageCount << std::endl;

        // write all elements linearly
        for(unsigned page_i = 0; page_i < pageCount; ++page_i)
        {
            for(unsigned j = 0; j < elementsPerPage; ++j)
            {
                file << std::fixed << std::setprecision(20) << data.at(j + page_i * elementsPerPage) << " , ";
            }
        }
        file << std::endl;
    }

    void writeMatlabOutput(
        std::filesystem::path const experimentPath,
        std::vector<float> const ase,
        std::vector<unsigned> const N_rays,
        std::vector<double> const expectedValues,
        unsigned const numberOfSamples,
        unsigned const numberOfLevels)
    {
        std::ofstream aseFile;
        std::ofstream raysFile;
        std::ofstream expectedValuesFile;
        unsigned const samplesPerLevel = numberOfSamples / numberOfLevels;

        aseFile.open(experimentPath / "phi_ASE.txt");
        write3dMatrix(ase, aseFile, samplesPerLevel, numberOfLevels, 1);
        aseFile.close();

        raysFile.open(experimentPath / "N_rays.txt");
        write3dMatrix(N_rays, raysFile, samplesPerLevel, numberOfLevels, 1);
        raysFile.close();

        expectedValuesFile.open(experimentPath / "mse_values.txt");
        write3dMatrix(expectedValues, expectedValuesFile, samplesPerLevel, numberOfLevels, 1);
        expectedValuesFile.close();
    }

} // namespace hase::utils
