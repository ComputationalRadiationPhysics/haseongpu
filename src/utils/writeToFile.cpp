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


#include <utils/writeToFile.hpp>

#include <ctime> /* time */
#include <fstream>
#include <iomanip> /* setw(), setfill */
#include <sstream> /* std::stringstream */
#include <string> /* std::to_string */
#include <vector>

namespace hase::utils
{
    namespace fs = std::filesystem;

    int writeValueToFile(
        float const value,
        fs::path const path,
        std::string const indexName1,
        int const index1,
        std::string const indexName2,
        int const index2)
    {
        std::stringstream filenameStream;
        filenameStream << path << indexName1 << "_" << std::setfill('0') << std::setw(3) << index1 << "_" << indexName2
                       << std::setfill('0') << std::setw(6) << index2;

        std::ofstream oFile;
        oFile.open(filenameStream.str());
        oFile << value << std::endl;
        oFile.close();

        return 0;
    }

    void writeVectorToFile(std::vector<double> v, fs::path filename)
    {
        // Add time to filename
        time_t currentTime;
        time(&currentTime);
        filename += "_";
        filename += std::to_string((int) currentTime);
        filename += ".dat";

        // Init filestream
        std::ofstream file;
        file.open(filename);

        // Write vector data
        for(std::vector<double>::iterator it = v.begin(); it != v.end(); ++it)
        {
            file << *it << std::endl;
        }
        file.close();
    }

} // namespace hase::utils
