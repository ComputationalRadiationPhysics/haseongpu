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
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @author Marius Melzer
 * @licence GPLv3
 *
 */

#pragma once
#include <core/logging.hpp>
#include <core/mesh.hpp>
#include <core/types.hpp>
#include <parse/CmdOptionsMap.hpp>

#include <fstream>
#include <string> /* string */
#include <vector>

namespace hase::parse
{

    /**
     * @brief Parses a given file(filename) line by line.
     *        Each line should contain just one value
     *        and the value should be a number (short, unsigned,
     *        int, float, double).
     *
     * @param filename file to parse
     * @return vector that contains the parsed values
     **/
    template<class T>
    std::vector<T> fileToVector(fs::path filename)
    {
        std::ifstream fileStream;

        fileStream.open(filename);

        if(fileStream.is_open())
        {
            std::vector<T> v;
            T value = 0.0;
            while(fileStream.good())
            {
                fileStream >> value;
                if(alpaka::math::isnan(value))
                {
                    hase::core::dout(V_ERROR) << "NAN in input data: " << filename << std::endl;
                    exit(1);
                }
                v.push_back(value);
            }
            v.pop_back();
            fileStream.close();
            return v;
        }
        else
        {
            hase::core::dout(V_ERROR) << "Can't open file " << filename << std::endl;
            exit(1);
        }
    }

    /**
     * @brief Parses just one line(value) of a given file(filename).
     *        The value should be a number (short, unsigned,
     *        int, float, double).
     *
     * @param filename file to parse
     * @return the value that was parsed
     **/
    template<class T>
    T fileToValue(fs::path filename)
    {
        std::ifstream fileStream;

        fileStream.open(filename);
        if(fileStream.is_open())
        {
            T value;
            fileStream >> value;
            fileStream.close();
            return value;
        }
        else
        {
            hase::core::dout(V_ERROR) << "Can't open file " << filename << std::endl;
            exit(1);
        }
    }

    hase::parse::CmdOptionsMap parseCommandLine(int const argc, char** argv);


    void printCommandLine(hase::parse::CmdOptionsMap const&);

    template<alpaka::onHost::concepts::Device T_Device>
    hase::core::DeviceMeshContainer<T_Device> createMesh(
        T_Device& device,
        std::vector<unsigned> const& triangleIndices,
        unsigned const numberOfTriangles,
        unsigned const numberOfLevels,
        unsigned const numberOfPoints,
        float const thicknessOfPrism,
        std::vector<double>& pointsVector,
        std::vector<double>& xOfTriangleCenter,
        std::vector<double>& yOfTriangleCenter,
        std::vector<unsigned>& positionsOfNormalVectors,
        std::vector<double>& xOfNormals,
        std::vector<double>& yOfNormals,
        std::vector<int>& forbiddenVector,
        std::vector<int>& neighborsVector,
        std::vector<float>& surfacesVector,
        std::vector<double>& betaValuesVector,
        std::vector<double>& betaCells,
        std::vector<unsigned>& cellTypes,
        std::vector<float>& refractiveIndices,
        std::vector<float>& reflectivities,
        std::vector<float>& totalReflectionAngles,
        float const nTot,
        float const crystalFluorescence,
        unsigned const cladNumber,
        double const cladAbsorption);
    hase::core::HostMesh createHostMeshFromFile(fs::path rootPath);
    int constructSimulationContextFromOptionsMap(
        hase::parse::CmdOptionsMap& vm,
        hase::core::ExperimentParameters& experiment,
        hase::core::ComputeParameters& compute,
        std::optional<hase::core::HostMesh>& hostMesh,
        hase::core::Result& result);
    hase::parse::CmdOptionsMap constructOptionsMapFromFile(std::filesystem::path const& cfgPath);
    int pythonParse(
        hase::core::ExperimentParameters& experiment,
        hase::core::ComputeParameters& compute,
        hase::core::HostMesh& host_mesh,
        hase::core::Result& result);

    int parse(
        int const argc,
        char** argv,
        hase::core::ExperimentParameters& experiment,
        hase::core::ComputeParameters& compute,
        std::optional<hase::core::HostMesh>& hostMesh,
        hase::core::Result& result);

} // namespace hase::parse
