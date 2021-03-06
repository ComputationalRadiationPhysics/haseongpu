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
 * @author Marius Melzer
 * @licence GPLv3
 *
 */

#pragma once
#include <string>  /* string */
#include <iostream>
#include <fstream> /* ifstream */
#include <cstdlib> /* atof */
#include <vector> 

#include <logging.hpp>
#include <mesh.hpp>
#include <nan_fix.hpp>

enum DeviceMode { NO_DEVICE_MODE, GPU_DEVICE_MODE, CPU_DEVICE_MODE};
enum ParallelMode { NO_PARALLEL_MODE, THREADED_PARALLEL_MODE, MPI_PARALLEL_MODE };


/**
 * @brief Parses a given file(filename) line by line.
 *        Each line should contain just one value
 *        and the value should be a number (short, unsigned,
 *        int, float, double).
 *
 * @param filename file to parse
 * @param vector contains the parsed values 
 *
 * @return 1 if parsing was succesful (file can be opened)
 * @return 0 otherwise
 **/
template <class T>
int fileToVector(const std::string filename, std::vector<T> *v){
  std::string line;
  std::ifstream fileStream; 
  T value = 0.0;

  fileStream.open(filename.c_str());
  if(fileStream.is_open()){
    while(fileStream.good()){
      std::getline(fileStream, line);
      value = (T) atof(line.c_str());
      if(isNaN(value)){
	dout(V_ERROR) << "NAN in input data: " << filename << std::endl;
	exit(1);
      }
      v->push_back(value);
    }

  }
  else{
    dout(V_ERROR) << "Can't open file " << filename << std::endl;
    fileStream.close();
    return 1;
  }
  v->pop_back();
  fileStream.close();
  return 0;
  
}
/**
 * @brief Parses just one line(value) of a given file(filename).
 *        The value should be a number (short, unsigned,
 *        int, float, double).
 *
 * @param filename file to parse
 * @param value is the value which was parsed 
 *
 * @return 1 if parsing was succesful (file can be opened)
 * @return 0 otherwise
 **/
template <class T>
int fileToValue(const std::string filename, T &value){
  std::string line;
  std::ifstream fileStream; 

  fileStream.open(filename.c_str());
  if(fileStream.is_open()){
      std::getline(fileStream, line);
      value = (T) atof(line.c_str());
  }
  else{
    dout(V_ERROR) << "Can't open file " << filename << std::endl;
    fileStream.close();
    return 1;
  }
  fileStream.close();
  return 0;
  
}


void parseCommandLine(
    const int argc,
    char** argv,
    unsigned *raysPerSample,
    unsigned *maxRaysPerSample,
    std::string *inputPath,
    bool *writeVtk,
    DeviceMode *deviceMode,
    ParallelMode *parallelMode,
    bool *useReflections,
    unsigned *maxgpus,
    int *minSample_i,
    int *maxSample_i,
    unsigned *maxRepetitions,
    std::string *outputPath,
    double *mseThreshold,
    unsigned *lambdaResolution
    );

void printCommandLine(
    unsigned raysPerSample,
    unsigned maxRaysPerSample,
    std::string inputPath,
    bool writeVtk,
    std::string compareLocation,
    const DeviceMode deviceMode,
    const ParallelMode parallelMode,
    bool useReflections,
    unsigned maxgpus,
    int minSample_i,
    int maxSample_i,
    unsigned maxRepetitions,
    std::string outputPath,
    double mseThreshold
    );

int checkParameterValidity(
    const int argc,
    const unsigned raysPerSample,
    unsigned *maxRaysPerSample,
    const std::string inputPath,
    const unsigned deviceCount,
    const DeviceMode deviceMode,
    const ParallelMode parallelMode,
    unsigned *maxgpus,
    const int minSample_i,
    const int maxSample_i,
    const unsigned maxRepetitions,
    const std::string outputPath,
    double *mseThreshold
    );

void checkSampleRange(
	int* minSampleRange,
	int* maxSampleRange,
	const unsigned numberOfSamples
	);

std::vector<Mesh> parseMesh(std::string rootPath,
			    std::vector<unsigned> devices,
			    unsigned maxGpus);


