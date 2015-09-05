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

// STL
#include <string>  /* string */
#include <vector>

// Boost
#include <boost/filesystem/fstream.hpp> /* fs::fstream */
#include <boost/filesystem/path.hpp> /* fs::path */
#include <boost/program_options/variables_map.hpp>

// HASEonGPU
#include <logging.hpp>
#include <mesh.hpp>
#include <nan_fix.hpp>
#include <types.hpp>


namespace fs = boost::filesystem;
namespace po = boost::program_options;

typedef std::map<std::string, po::variable_value> Modifiable_variables_map;


/**
 * @brief Parses a given file(filename) line by line.
 *        Each line should contain just one value
 *        and the value should be a number (short, unsigned,
 *        int, float, double).
 *
 * @param filename file to parse
 * @return vector that contains the parsed values
 **/
template <class T>
std::vector<T> fileToVector(fs::path filename){
    fs::ifstream fileStream;

    fileStream.open(filename);

    if(fileStream.is_open()){
        std::vector<T> v;
        T value = 0.0;
        while(fileStream.good()){
            fileStream >> value;
            if(isNaN(value)){
                dout(V_ERROR) << "NAN in input data: " << filename << std::endl;
                exit(1);
            }
            v.push_back(value);
        }
        v.pop_back();
        fileStream.close();
        return v;

    }
    else{
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
template <class T>
T fileToValue(fs::path filename){
    fs::ifstream fileStream;

    fileStream.open(filename);
    if(fileStream.is_open()){
        T value;
        fileStream >> value;
        fileStream.close();
        return value;
    }
    else{
        dout(V_ERROR) << "Can't open file " << filename << std::endl;
        exit(1);
    }
}


po::variables_map parseCommandLine(const int argc, char** argv);


void printCommandLine(const Modifiable_variables_map);


Modifiable_variables_map checkParameterValidity(Modifiable_variables_map, unsigned);


Modifiable_variables_map checkSampleRange(
    Modifiable_variables_map vm,
    const unsigned numberOfSamples
    );


// std::vector<Mesh> parseMesh(const fs::path rootPath,
// 			    std::vector<unsigned> devices);


int parse( const int argc,
	   char** argv,
	   ExperimentParameters& experiment,
	   ComputeParameters& compute,
	   std::vector<Mesh>& mesh,
	   Result& result);
