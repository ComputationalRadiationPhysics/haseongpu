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



#include <fstream> /* ofstream */
#include <vector> /* vector */
#include <iomanip> /* std::setprecision() */
#include <string>
#include <time.h> /* time, time_t */
#include <sstream> /* std::stringstream */

#include <write_to_vtk.hpp>
#include <logging.hpp>
#include <mesh.hpp>


/** 
 * @brief takes data and creates a nice VTK-file usable with paraview
 *
 * @param 
 */
int writeToVtk(const Mesh& mesh,
	       const std::vector<double> data,
	       const std::string pfilename,
	       const unsigned raysPerSample,
	       const unsigned maxRaysPerSample,
	       const float expectationThreshold,
	       const bool useReflections,
	       const float runtime,
         const std::string vtkType){

  const double* vertexCoordinates = mesh.points.toArray();
  const unsigned* triangles       = mesh.trianglePointIndices.toArray();
  const float    thicknessOfLevel = mesh.thickness;
  const unsigned verticesPerLevel = mesh.numberOfPoints;
  const unsigned trianglesPerLevel= mesh.numberOfTriangles;
  const unsigned numberOfLevels   = mesh.numberOfLevels;
  const unsigned numberOfCells    = trianglesPerLevel*(numberOfLevels-1);
  const unsigned numberOfVertices = numberOfLevels*verticesPerLevel;

  // Construct experiment information
  unsigned r = useReflections ? mesh.getMaxReflections() : 0;
  
  std::stringstream experimentStream;
  experimentStream << "RAYS=" << raysPerSample << " MAXRAYS=" << maxRaysPerSample << " REFLECTIONS=" << r << " EXPECTATION=" << expectationThreshold << " RUNTIME=" << runtime;

  // Add time to filename
  time_t currentTime;
  time(&currentTime);
  std::stringstream filenameStream;
  filenameStream  << pfilename << "_" << (int) currentTime << ".vtk";

  dout(V_INFO) << "Write experiment data to vtk-file " << filenameStream.str() << std::endl;

  std::ofstream vtkFile;
  vtkFile.open(filenameStream.str().c_str());

  // Write header of vtk file
  vtkFile << "# vtk DataFile Version 2.0" << std::endl;
  vtkFile << experimentStream.str() << std::endl;
  vtkFile << "ASCII" << std::endl;

  // Write point data
  vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
  vtkFile << "POINTS " << verticesPerLevel*numberOfLevels <<  " float" << std::endl;
  for(unsigned level_i=0; level_i < numberOfLevels ; ++level_i){
    for(unsigned point_i=0; point_i < verticesPerLevel; ++point_i){
      vtkFile << std::fixed << std::setprecision(6) << vertexCoordinates[point_i] << " " << vertexCoordinates[point_i + verticesPerLevel] << " " << level_i * thicknessOfLevel << std::endl;
    }
  }

  // Write cell data
  vtkFile << "CELLS" << " " << numberOfCells << " " << numberOfCells * 7 << std::endl;
  for(unsigned level_i=0; level_i < (numberOfLevels - 1); ++level_i){
    for(unsigned triangle_i=0; triangle_i < trianglesPerLevel; ++triangle_i){
      vtkFile << "6 " 
        << level_i     * verticesPerLevel + triangles[triangle_i] << " "
        << level_i     * verticesPerLevel + triangles[trianglesPerLevel + triangle_i] << " "
        << level_i     * verticesPerLevel + triangles[2 * trianglesPerLevel + triangle_i] << " "
        << (level_i+1) * verticesPerLevel + triangles[triangle_i] << " "
        << (level_i+1) * verticesPerLevel + triangles[trianglesPerLevel + triangle_i] << " "
        << (level_i+1) * verticesPerLevel + triangles[2 * trianglesPerLevel + triangle_i] << std::endl;

    }

  }

  // Write cell type
  vtkFile << "CELL_TYPES " << numberOfCells << std::endl;
  for(unsigned i=0; i < numberOfCells; ++i){
    // 13 is the VTK type for this kind of cell (prism/wedge)
    vtkFile << "13" << std::endl;
  }

  // Write data
  vtkFile << vtkType << " " << numberOfVertices << std::endl;
  vtkFile << "SCALARS scalars float 1" << std::endl;
  vtkFile << "LOOKUP_TABLE default" << std::endl;
  for(unsigned i=0; i < numberOfVertices; ++i){
    vtkFile << std::fixed << std::setprecision(6) << data.at(i) << std::endl;
  }

  vtkFile.close();

  return 0;
}


int writePrismToVtk(
    const Mesh& mesh,
    const std::vector<double> prismData,
    const std::string pfilename,
    const unsigned raysPerSample,
    const unsigned maxRaysPerSample,
    const float expectationThreshold,
    const bool useReflections,
    const float runtime){


  return writeToVtk(
      mesh,
      prismData,
      pfilename,
      raysPerSample,
      maxRaysPerSample,
      expectationThreshold,
      useReflections,
      runtime,
      "CELL_DATA");

}

int writePointsToVtk(
    const Mesh& mesh,
    const std::vector<double> prismData,
    const std::string pfilename,
    const unsigned raysPerSample,
    const unsigned maxRaysPerSample,
    const float expectationThreshold,
    const bool useReflections,
    const float runtime){


 return writeToVtk(
      mesh,
      prismData,
      pfilename,
      raysPerSample,
      maxRaysPerSample,
      expectationThreshold,
      useReflections,
      runtime,
      "POINT_DATA");
}



std::vector<double> compareVtk(std::vector<double> compare, std::string filename){
  std::ifstream filestream;
  std::string line;
  bool foundLine = false;
  double value = 0;
  double diff = 0;
  unsigned ase_i = 0;
  double minDiff = 10000; // should be enough
  double maxDiff = 0;
  double totalDiff = 0;
  double aseTotal = 0;
  double smallDiff = 10;

  // No compare vtk was given
  if(!filename.compare("")){
    return std::vector<double>();
  }
  dout(V_INFO) << "Compare solution with " << filename << std::endl;

  for(unsigned i = 0; i < compare.size(); ++i){
    aseTotal += compare.at(i);
  }

  filestream.open(filename.c_str(), std::ifstream::in);

  if(filestream.is_open()){
    while(filestream.good()){
      std::getline(filestream, line);
      std::size_t found = line.find("LOOKUP_TABLE default");
      if(found != std::string::npos){
        foundLine = true;
        std::getline(filestream, line);
      }
      if(foundLine){
        if(ase_i == compare.size())
          break;
        value = (double) atof(line.c_str());

        if(abs(value) > abs(compare.at(ase_i)))
          diff = (abs(value / compare.at(ase_i)) - 1) * 100;
        else
          diff = (abs(compare.at(ase_i) / value) - 1) * 100;

        totalDiff += diff;

        if(diff >= maxDiff)
          maxDiff = diff;

        if(diff <= minDiff)
          minDiff = diff;

        if(diff >= smallDiff){
          dout(V_WARNING) << "ASE relative difference[" << ase_i << "]: " << diff  << "%" << "[" << compare.at(ase_i) << ", " <<     value  << "]"<<" > " << smallDiff << "%" << std::endl;
        }
        compare.at(ase_i) = diff;
        ase_i++;

      }

    }

  }
  else{
    dout(V_WARNING) << "Can't open file " << filename << " for comparison" << std::endl;
    return std::vector<double>();
  }

  dout(V_STAT) << "ASE max. difference: " << maxDiff << "%" << std::endl;
  dout(V_STAT) << "ASE min. difference: " << minDiff << "%" << std::endl;
  dout(V_STAT) << "ASE tot. avg difference: " << totalDiff / compare.size() << "%" << std::endl;
  filestream.close();
  return compare;
}
