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


#include <core/logging.hpp>
#include <core/mesh.hpp>
#include <utils/writeToVtk.hpp>

#include <filesystem>
#include <fstream>
#include <iomanip> /* std::fixed, std::setprecision() */
#include <string> /* std::to_string() */
#include <vector> /* vector */

namespace hase::utils
{

    /**
     * @brief takes data and creates a nice VTK-file usable with paraview
     *
     * @param
     */
    int writeToVtk(
        hase::core::HostMesh const& mesh,
        std::vector<double> const data,
        std::filesystem::path filename,
        unsigned const raysPerSample,
        unsigned const maxRaysPerSample,
        float const expectationThreshold,
        bool const useReflections,
        float const runtime,
        std::string const vtkType)
    {
        hase::core::dout(V_INFO) << "Write experiment data to vtk-file " << filename << std::endl;

        std::ofstream vtkFile;
        vtkFile.open(filename);

        vtkFile << "# vtk DataFile Version 2.0" << std::endl;
        vtkFile << "RAYS=" << raysPerSample << " MAXRAYS=" << maxRaysPerSample << " REFLECTIONS=0"
                << " EXPECTATION=" << expectationThreshold << " RUNTIME=" << runtime
                << " REQUESTED_REFLECTIONS=" << useReflections << std::endl;
        vtkFile << "ASCII" << std::endl;
        vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

        vtkFile << "POINTS " << mesh.numberOfPoints << " float" << std::endl;
        for(unsigned point = 0u; point < mesh.numberOfPoints; ++point)
        {
            vtkFile << std::fixed << std::setprecision(6) << mesh.points.at(point) << " "
                    << mesh.points.at(point + mesh.numberOfPoints) << " "
                    << mesh.points.at(point + 2u * mesh.numberOfPoints) << std::endl;
        }

        vtkFile << "CELLS " << mesh.numberOfCells << " " << mesh.numberOfCells * 7u << std::endl;
        for(unsigned cell = 0u; cell < mesh.numberOfCells; ++cell)
        {
            vtkFile << "6";
            for(unsigned localVertex = 0u; localVertex < hase::core::prism6VertexCount; ++localVertex)
            {
                vtkFile << " " << mesh.cellPointIndices.at(cell * hase::core::prism6VertexCount + localVertex);
            }
            vtkFile << std::endl;
        }

        vtkFile << "CELL_TYPES " << mesh.numberOfCells << std::endl;
        for(unsigned cell = 0u; cell < mesh.numberOfCells; ++cell)
        {
            vtkFile << mesh.cellTypes.at(cell) << std::endl;
        }

        vtkFile << vtkType << " " << data.size() << std::endl;
        vtkFile << "SCALARS scalars float 1" << std::endl;
        vtkFile << "LOOKUP_TABLE default" << std::endl;
        for(double const value : data)
        {
            vtkFile << std::fixed << std::setprecision(6) << value << std::endl;
        }

        vtkFile.close();
        return 0;
    }

    int writePrismToVtk(
        hase::core::HostMesh const& mesh,
        std::vector<double> const prismData,
        std::filesystem::path const filename,
        unsigned const raysPerSample,
        unsigned const maxRaysPerSample,
        float const expectationThreshold,
        bool const useReflections,
        float const runtime)
    {
        return writeToVtk(
            mesh,
            prismData,
            filename,
            raysPerSample,
            maxRaysPerSample,
            expectationThreshold,
            useReflections,
            runtime,
            "CELL_DATA");
    }

    int writePointsToVtk(
        hase::core::HostMesh const& hostMesh,
        std::vector<double> const prismData,
        std::filesystem::path const filename,
        unsigned const raysPerSample,
        unsigned const maxRaysPerSample,
        float const expectationThreshold,
        bool const useReflections,
        float const runtime)
    {
        return writeToVtk(
            hostMesh,
            prismData,
            filename,
            raysPerSample,
            maxRaysPerSample,
            expectationThreshold,
            useReflections,
            runtime,
            "POINT_DATA");
    }

    std::vector<double> compareVtk(std::vector<double> compare, std::filesystem::path const filename)
    {
        std::ifstream filestream;
        std::string line;
        bool foundLine = false;
        double value = 0;
        double diff = 0;
        unsigned ase_i = 0;
        double minDiff = 10000; // should be enough
        double maxDiff = 0;
        double totalDiff = 0;
        double smallDiff = 10;

        // No compare vtk was given
        if(!filename.compare(std::filesystem::path("")))
        {
            return std::vector<double>();
        }
        hase::core::dout(V_INFO) << "Compare solution with " << filename << std::endl;

        filestream.open(filename, std::ifstream::in);

        if(filestream.is_open())
        {
            while(filestream.good())
            {
                std::getline(filestream, line);
                std::size_t found = line.find("LOOKUP_TABLE default");
                if(found != std::string::npos)
                {
                    foundLine = true;
                    std::getline(filestream, line);
                }
                if(foundLine)
                {
                    if(ase_i == compare.size())
                        break;
                    value = (double) atof(line.c_str());

                    if(alpaka::math::abs(value) > alpaka::math::abs(compare.at(ase_i)))
                        diff = (alpaka::math::abs(value / compare.at(ase_i)) - 1) * 100;
                    else
                        diff = (alpaka::math::abs(compare.at(ase_i) / value) - 1) * 100;

                    totalDiff += diff;

                    if(diff >= maxDiff)
                        maxDiff = diff;

                    if(diff <= minDiff)
                        minDiff = diff;

                    if(diff >= smallDiff)
                    {
                        hase::core::dout(V_WARNING)
                            << "ASE relative difference[" << ase_i << "]: " << diff << "%" << "[" << compare.at(ase_i)
                            << ", " << value << "]" << " > " << smallDiff << "%" << std::endl;
                    }
                    compare.at(ase_i) = diff;
                    ase_i++;
                }
            }
        }
        else
        {
            hase::core::dout(V_WARNING) << "Can't open file " << filename << " for comparison" << std::endl;
            return std::vector<double>();
        }

        hase::core::dout(V_STAT) << "ASE max. difference: " << maxDiff << "%" << std::endl;
        hase::core::dout(V_STAT) << "ASE min. difference: " << minDiff << "%" << std::endl;
        hase::core::dout(V_STAT) << "ASE tot. avg difference: " << totalDiff / compare.size() << "%" << std::endl;
        filestream.close();
        return compare;
    }

} // namespace hase::utils
