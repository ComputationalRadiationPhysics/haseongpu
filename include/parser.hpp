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


/**
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @author Marius Melzer
 * @licence GPLv3
 *
 */

#pragma once
#include <logging.hpp>
#include <mesh.hpp>
#include <nan_fix.hpp>
#include <types.hpp>

#include <string> /* string */
#include <vector>
#include <fstream>
#include <CmdOptionsMap.hpp>

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
            if(isNaN(value))
            {
                dout(V_ERROR) << "NAN in input data: " << filename << std::endl;
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
        dout(V_ERROR) << "Can't open file " << filename << std::endl;
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
        dout(V_ERROR) << "Can't open file " << filename << std::endl;
        exit(1);
    }
}

CmdOptionsMap parseCommandLine(int const argc, char** argv);


void printCommandLine(CmdOptionsMap const &);

Mesh createMesh(
    std::vector<unsigned> const& triangleIndices,
    unsigned numberOfTriangles,
    unsigned numberOfLevels,
    unsigned numberOfPoints,
    float thicknessOfPrism,
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
    float nTot,
    float crystalFluorescence,
    unsigned cladNumber,
    double cladAbsorption);


std::vector<Mesh> parseMesh(fs::path const rootPath, std::vector<unsigned> devices);

int pythonParse(
    ExperimentParameters& experiment,
    ComputeParameters& compute,
    HostMesh& host_mesh,
    std::vector<Mesh>& mesh,
    Result& result);

int parse(
    int const argc,
    char** argv,
    ExperimentParameters& experiment,
    ComputeParameters& compute,
    std::vector<Mesh>& mesh,
    Result& result);
