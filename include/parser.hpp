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
#include <nan_fix.hpp>
#include <types.hpp>
#include <mesh.hpp>     


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

int parse( const int argc,
	   char** argv,
	   ExperimentParameters& experiment,
	   ComputeParameters& compute,
	   Result& result);




template <class T, class B, class E>
void assertRange(const std::vector<T> &v, const B minElement,const E maxElement, const bool equals){
  if(equals){
    assert(*std::min_element(v.begin(),v.end()) == minElement);
    assert(*std::max_element(v.begin(),v.end()) == maxElement);
  }else{
    assert(*std::min_element(v.begin(),v.end()) >= minElement);
    assert(*std::max_element(v.begin(),v.end()) <= maxElement);
  }
}


template <class T, class B>
void assertMin(const std::vector<T>& v,const  B minElement,const bool equals){
  if(equals){
    assert(*std::min_element(v.begin(),v.end()) == minElement);
  }else{
    assert(*std::min_element(v.begin(),v.end()) >= minElement);
  }
}


/**
 * @brief fills a device mesh with the correct datastructures
 *
 */
template <typename T_Acc, typename T_Dev>
Mesh<T_Acc, T_Dev> createMesh(const std::vector<unsigned> &triangleIndices, 
			      const unsigned numberOfTriangles, 
			      const unsigned numberOfLevels,
			      const unsigned numberOfPoints, 
			      const float thicknessOfPrism,
			      std::vector<double> &pointsVector, 
			      std::vector<double> &xOfTriangleCenter, 
			      std::vector<double> &yOfTriangleCenter, 
			      std::vector<unsigned> &positionsOfNormalVectors,
			      std::vector<double> &xOfNormals, 
			      std::vector<double> &yOfNormals,
			      std::vector<int> &forbiddenVector, 
			      std::vector<int> &neighborsVector, 
			      std::vector<float> &surfacesVector,
			      std::vector<double> &betaValuesVector,
			      std::vector<double> &betaCells,
			      std::vector<unsigned> &cellTypes,
			      std::vector<float> & refractiveIndices,
			      std::vector<float> & reflectivities,
			      const float nTot,
			      const float crystalFluorescence,
			      const unsigned cladNumber,
			      const double cladAbsorption,
			      T_Dev &dev) {

    // GPU variables
    double totalSurface = 0.;

    for(unsigned i=0;i<numberOfTriangles;++i){
	totalSurface+=double(surfacesVector.at(i));	
    }

    // Vector Preprocessing
    std::vector<double> hostNormalVec(xOfNormals.begin(), xOfNormals.end());
    hostNormalVec.insert(hostNormalVec.end(),yOfNormals.begin(),yOfNormals.end());
    std::vector<double> hostCenters(xOfTriangleCenter.begin(), xOfTriangleCenter.end());
    hostCenters.insert(hostCenters.end(),yOfTriangleCenter.begin(),yOfTriangleCenter.end());
    std::vector<float> totalReflectionAngles(refractiveIndices.size()/2,0);
    for(unsigned i=0;i<refractiveIndices.size();i+=2){
	totalReflectionAngles.at(i/2) = (180. / M_PI *  asin(refractiveIndices.at(i+1) / refractiveIndices.at(i)));
    }

    return Mesh<T_Acc, T_Dev> ( cladAbsorption,
				totalSurface,
				thicknessOfPrism,
				nTot,
				crystalFluorescence,
				numberOfTriangles,
				numberOfLevels,
				numberOfTriangles * (numberOfLevels-1),
				numberOfPoints,
				numberOfPoints * numberOfLevels,
				cladNumber,
				pointsVector,
				hostNormalVec,
				betaValuesVector,
				hostCenters,
				surfacesVector,
				forbiddenVector,
				betaCells,
				cellTypes,
				refractiveIndices,
				reflectivities,
				totalReflectionAngles,
				triangleIndices,
				neighborsVector,
				positionsOfNormalVectors,
				dev);

}


template <typename T_Acc, typename T_Dev>
Mesh<T_Acc, T_Dev> parseMesh( fs::path const rootPath,
			      T_Dev &dev) {

    // Parse experimentdata from files
    std::vector<unsigned> triangleNormalPoint  = fileToVector<unsigned>(rootPath / "triangleNormalPoint.txt");
    std::vector<double> betaVolume             = fileToVector<double>(rootPath / "betaVolume.txt");
    std::vector<int> forbiddenEdge             = fileToVector<int>(rootPath / "forbiddenEdge.txt");
    std::vector<int> triangleNeighbors         = fileToVector<int>(rootPath / "triangleNeighbors.txt");
    std::vector<double> triangleNormalsX       = fileToVector<double>(rootPath / "triangleNormalsX.txt");
    std::vector<double> triangleNormalsY       = fileToVector<double>(rootPath / "triangleNormalsY.txt");
    std::vector<double> triangleCenterX        = fileToVector<double>(rootPath / "triangleCenterX.txt");
    std::vector<double> triangleCenterY        = fileToVector<double>(rootPath / "triangleCenterY.txt");
    std::vector<double> points                 = fileToVector<double>(rootPath / "points.txt");
    std::vector<unsigned> trianglePointIndices = fileToVector<unsigned>(rootPath / "trianglePointIndices.txt");
    std::vector<float>  triangleSurfaces       = fileToVector<float>(rootPath / "triangleSurfaces.txt");
    unsigned numberOfPoints                    = fileToValue<unsigned>(rootPath / "numberOfPoints.txt");
    unsigned numberOfTriangles                 = fileToValue<unsigned>(rootPath / "numberOfTriangles.txt");
    unsigned numberOfLevels                    = fileToValue<unsigned>(rootPath / "numberOfLevels.txt");
    float thickness                            = fileToValue<float>(rootPath / "thickness.txt");
    float nTot                                 = fileToValue<float>(rootPath / "nTot.txt");
    float crystalTFluo                         = fileToValue<float>(rootPath / "crystalTFluo.txt");
    unsigned claddingNumber                    = fileToValue<unsigned>(rootPath / "claddingNumber.txt");
    double claddingAbsorption                  = fileToValue<double>(rootPath / "claddingAbsorption.txt");
    std::vector<double>  betaCells             = fileToVector<double>(rootPath / "betaCells.txt");
    std::vector<unsigned>  claddingCellTypes   = fileToVector<unsigned>(rootPath / "claddingCellTypes.txt");
    std::vector<float>  refractiveIndices      = fileToVector<float>(rootPath / "refractiveIndices.txt");
    std::vector<float>  reflectivities         = fileToVector<float>(rootPath / "reflectivities.txt");

    // assert input sizes
    assert(numberOfPoints == (points.size() / 2));
    assert(numberOfTriangles == trianglePointIndices.size() / 3);
    assert(triangleNormalPoint.size() == numberOfTriangles * 3);
    assert(triangleCenterY.size() == numberOfTriangles);
    assert(triangleCenterX.size() == numberOfTriangles);
    assert(triangleSurfaces.size() == numberOfTriangles);
    assert(betaVolume.size() == numberOfTriangles * (numberOfLevels-1));
    assert(triangleNormalsX.size() == numberOfTriangles * 3);
    assert(triangleNormalsY.size() == numberOfTriangles * 3);
    assert(trianglePointIndices.size() == numberOfTriangles * 3);
    assert(forbiddenEdge.size() == numberOfTriangles * 3);
    assert(triangleNeighbors.size() == numberOfTriangles * 3);
    assert(betaCells.size() == numberOfPoints * numberOfLevels);
    assert(claddingCellTypes.size()== numberOfTriangles);
    assert(refractiveIndices.size() == 4);
    assert(reflectivities.size() == (refractiveIndices.size()/2) * numberOfTriangles);
    assert(claddingCellTypes.size() == numberOfTriangles);

    // assert input data validity
    assertRange(triangleNormalPoint,0u,unsigned(numberOfPoints-1),true);
    assertMin(betaVolume,0,false);
    assertRange(forbiddenEdge,-1,2,true);
    assertRange(triangleNeighbors,-1,int(numberOfTriangles-1),true);
    assertRange(triangleNormalsX,-1,1,false);
    assertRange(triangleNormalsY,-1,1,false);
    assertRange(triangleCenterX,*std::min_element(points.begin(),points.end()),*std::max_element(points.begin(),points.end()),false);
    assertRange(triangleCenterY,*std::min_element(points.begin(),points.end()),*std::max_element(points.begin(),points.end()),false);
    assertRange(trianglePointIndices,0u,unsigned(numberOfPoints-1),true);
    assertMin(triangleSurfaces,0,false);
    assertMin(betaCells,0,false);
    assertRange(refractiveIndices,0,5,false);
    assertRange(reflectivities,0,1,false);

    return createMesh<T_Acc, T_Dev>(trianglePointIndices,
				    numberOfTriangles,
				    numberOfLevels,
				    numberOfPoints,
				    thickness,
				    points,
				    triangleCenterX,
				    triangleCenterY,
				    triangleNormalPoint,
				    triangleNormalsX,
				    triangleNormalsY,
				    forbiddenEdge,
				    triangleNeighbors,
				    triangleSurfaces,
				    betaVolume,
				    betaCells,
				    claddingCellTypes,
				    refractiveIndices,
				    reflectivities,
				    nTot,
				    crystalTFluo,
				    claddingNumber,
				    claddingAbsorption,
				    dev);
}

