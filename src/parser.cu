/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 *
 * This file is part of HASENonGPU
 *
 * HASENonGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HASENonGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HASENonGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */



#include <string> /* string */
#include <vector> /* vector */
#include <algorithm>
#include <assert.h>

#include <types.hpp>
#include <logging.hpp> 
#include <mesh.hpp>
#include <parser.hpp>

void parseCommandLine(
    const int argc,
    char** argv,
    unsigned *raysPerSample,
    unsigned *maxRaysPerSample,
    std::string *inputPath,
    bool *writeVtk,
    std::string *compareLocation,
    RunMode *mode,
    bool *useReflections,
    unsigned *maxgpus,
    int *minSample_i,
    int *maxSample_i,
    unsigned *maxRepetitions,
    std::string *outputPath,
    double *mseThreshold,
    unsigned *maxLambdaResolution
    ) {

  std::vector<std::pair<std::string, std::string> > parameters;

  // Parse Commandline
  for (int i = 1; i < argc; ++i) {

    char* pos = strtok(argv[i], "=");
    std::pair < std::string, std::string > p(std::string(pos), std::string(""));
    pos = strtok(NULL, "=");
    if (pos != NULL) {
      p.second = std::string(pos);
    }
    parameters.push_back(p);
  }
  for (unsigned i = 0; i < parameters.size(); ++i) {
    std::pair < std::string, std::string > p = parameters.at(i);
    dout(V_INFO) << "arg[" << i << "]: (" << p.first << "," << p.second << ")" << std::endl;

    // Parse number of rays
    if (p.first == "--rays") {
      *raysPerSample = atoi(p.second.c_str());
    }

    if (p.first == "--maxrays"){
      *maxRaysPerSample = atoi(p.second.c_str());
    }

    if(p.first == "--input") {
      std::string temp_input(p.second);

      // Add slash at the end, if missing
      if ((temp_input)[temp_input.size() - 1] == 'w')
        temp_input.erase(temp_input.size() - 1, 1);
      else if (temp_input[temp_input.size() - 1] != '/')
        temp_input.append("/");

      *inputPath = temp_input;
    }

    if( p.first =="--output") {

      std::string temp_output(p.second);

      // Add slash at the end, if missing
      if ((temp_output)[temp_output.size() - 1] == 'w')
        temp_output.erase(temp_output.size() - 1, 1);
      else if (temp_output[temp_output.size() - 1] != '/')
        temp_output.append("/");

      *outputPath = temp_output;
    }

    if (p.first == "--write-vtk") {
      *writeVtk = true;
    }

    // Parse what vtk file to compare with
    if (p.first == "--compare") {
      *compareLocation = p.second;
    }

    if (p.first == "--runmode") {
      if (p.second == "threaded")
        *mode = GPU_THREADED;
      if (p.second == "cpu")
        *mode = CPU;
      if (p.second == "mpi")
        *mode = GPU_MPI;

    }

    if (p.first == "--reflection"){
      *useReflections = true;
    }

    if (p.first == "--maxgpus"){
      *maxgpus = atoi(p.second.c_str());
    }

    if (p.first == "--min_sample_i"){
      *minSample_i = atoi(p.second.c_str());
    }
    if (p.first == "--max_sample_i"){
      *maxSample_i = atoi(p.second.c_str());
    }

    if (p.first == "--verbosity"){
      verbosity = unsigned(atoi(p.second.c_str()));
    }

    if(p.first == "--repetitions"){
      *maxRepetitions = unsigned(atoi(p.second.c_str()));
    }

    if(p.first == "--mse-threshold"){
      *mseThreshold = float(atof(p.second.c_str()));
    }

    if(p.first == "--max-lambda-resolution"){
      *maxLambdaResolution = unsigned(atoi(p.second.c_str()));
    }


  }
}

int checkParameterValidity(
    const int argc,
    const unsigned raysPerSample,
    unsigned *maxRaysPerSample,
    const std::string inputPath,
    const unsigned deviceCount,
    const RunMode mode,
    unsigned *maxgpus,
    const int minSampleRange,
    const int maxSampleRange,
    const unsigned maxRepetitions,
    const std::string outputPath,
    double *mseThreshold
    ) {

  if (argc <= 1) {
    dout(V_ERROR) << "No commandline arguments found" << std::endl;
    dout(V_ERROR) << "Usage: ./calcPhiASE ARGS [OPTARGS]" << std::endl;
    dout(V_ERROR) << "" << std::endl;
    dout(V_ERROR) << "ARGS (required)" << std::endl;
    dout(V_ERROR) << "  --runmode=RUNMODE" << std::endl;
    dout(V_ERROR) << "  --rays=RAYS" << std::endl;
    dout(V_ERROR) << "  --maxrays=MAXRAYS" << std::endl;
    dout(V_ERROR) << "  --input=FOLDER" << std::endl;
    dout(V_ERROR) << "  --output=FOLDER" << std::endl;
    dout(V_ERROR) << "" << std::endl;
    dout(V_ERROR) << "OPTARGS (optional)" << std::endl;
    dout(V_ERROR) << "  --maxgpus=N" << std::endl;
    dout(V_ERROR) << "  --min_sample_i=<index of first sample>" << std::endl;
    dout(V_ERROR) << "  --max_sample_i=<index of last sample>" << std::endl;
    dout(V_ERROR) << "  --compare=<location of vtk-file to compare with>" << std::endl;
    dout(V_ERROR) << "  --verbosity=VERBOSITY_LEVEL" << std::endl;
    dout(V_ERROR) << "  --repetitions=MAX_REPETITIONS" << std::endl;
    dout(V_ERROR) << "  --mse_threshold=THRESHOLD" << std::endl;
    dout(V_ERROR) << "  --reflection" << std::endl;
    dout(V_ERROR) << "  --write-vtk" << std::endl;
    dout(V_ERROR) << "" << std::endl;
    dout(V_ERROR) << "Runmodes : cpu" << std::endl;
    dout(V_ERROR) << "           threaded" << std::endl;
    dout(V_ERROR) << "           mpi" << std::endl;
    dout(V_ERROR) << "Verbosity levels: 0 (quiet)" << std::endl; 
    dout(V_ERROR) << "                  1 (error)" << std::endl; 
    dout(V_ERROR) << "                  2 (warning)" << std::endl; 
    dout(V_ERROR) << "                  4 (info)" << std::endl; 
    dout(V_ERROR) << "                  8 (statistics)" << std::endl; 
    dout(V_ERROR) << "                 16 (debug)" << std::endl; 
    dout(V_ERROR) << "                 32 (progressbar)" << std::endl; 
    dout(V_ERROR) << "" << std::endl; 
    dout(V_ERROR) << "Please see README for more details!" << std::endl; 
    return 1;
  }
  if (mode == NONE) {
    dout(V_ERROR) << "Please specify the runmode with --runmode=MODE" << std::endl;
    return 1;
  }

  if (raysPerSample == 0) {
    dout(V_ERROR) << "Please specify the number of rays per sample Point with --rays=RAYS" << std::endl;
    return 1;
  }

  if (inputPath.size() == 0) {
    dout(V_ERROR) << "Please specify the experiment's location with --input=PATH_TO_EXPERIMENT" << std::endl;
    return 1;
  }

  if (outputPath.size() == 0) {
    dout(V_ERROR) << "Please specify the output location with --output=PATH_TO_EXPERIMENT" << std::endl;
    return 1;
  }


  if(*maxRaysPerSample < raysPerSample){
    dout(V_WARNING) << "maxRays < raysPerSample. Increasing maxRays to " << raysPerSample << " (will be non-adaptive!)" << std::endl;
    *maxRaysPerSample = raysPerSample;
  }

  if(*maxgpus > deviceCount){ dout(V_ERROR) << "You don't have so many devices, use --maxgpus=" << deviceCount << std::endl;
    return 1;
  }

  if(*maxgpus == 0){
    *maxgpus = deviceCount;
  }

  if(maxSampleRange < minSampleRange){
    dout(V_ERROR) << "maxSampleRange < minSampleRange!" << std::endl;
    return 1;
  }

  int samplesForNode = maxSampleRange-minSampleRange+1;
  if(samplesForNode < int(*maxgpus) && (minSampleRange != -1 || maxSampleRange != -1)){
    dout(V_WARNING) << "More GPUs requested than there are sample points. Number of used GPUs reduced to " << samplesForNode << std::endl;
     *maxgpus = unsigned(samplesForNode);
  }

  if(verbosity >= 64){
    verbosity = 63;
    dout(V_WARNING) << "Verbosity level should be between 0 (quiet) and 63 (all). Levels can be bitmasked together." << std::endl;
  }
  if(maxRepetitions < 1){
    dout(V_ERROR) << "At least 1 repetition is necessary!" << std::endl;
  }

  if(*mseThreshold == 0){
    *mseThreshold = 1000;
  }

  return 0;
}

void checkSampleRange(int* minSampleRange, int* maxSampleRange, const unsigned numberOfSamples){
  if(*minSampleRange == -1 && *maxSampleRange== -1){
    dout(V_WARNING) << "minSample_i/maxSample_i not set! Assuming a sample range of " << std::endl;
    dout(V_WARNING) << "0 to " << numberOfSamples-1 << std::endl;
    *minSampleRange = 0;
    *maxSampleRange = numberOfSamples-1;
    return;
  }

  if((*minSampleRange == -1 && *maxSampleRange != -1) || (*minSampleRange != -1 && *maxSampleRange == -1)){
    dout(V_ERROR) << "check minSample_i/maxSample_i! (Allowed Range from 0 to " << numberOfSamples << ")";
    exit(1);
  }

  if((*maxSampleRange >= int(numberOfSamples) || *maxSampleRange < (int)0)){
    dout(V_ERROR) << "maxSample_i is out of range! (There are only " << numberOfSamples << " samples)";
    exit(1);
  }

  if((*minSampleRange < -1 || *minSampleRange >= numberOfSamples)){
    dout(V_ERROR) << "minSample_i is out of range! (There are only " << numberOfSamples << " samples)";
    exit(1);
  }
}

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
 * See parseMultiGPU for details on the parameters
 */
Mesh createMesh(const std::vector<unsigned> &triangleIndices, 
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
		const double cladAbsorption
	       ) {

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

   Mesh mesh( cladAbsorption,
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
	     positionsOfNormalVectors);
  return mesh;

}

/**
 *
 */
std::vector<Mesh> parseMesh(std::string rootPath,
			    std::vector<unsigned> devices,
			    unsigned maxGpus) {

  std::vector<Mesh> meshs;

  // Experimentdata
  std::vector<double>  betaVolume;
  std::vector<double>  triangleNormalsX;
  std::vector<double>  triangleNormalsY;
  std::vector<unsigned>  trianglePointIndices;
  std::vector<int>  forbiddenEdge;
  std::vector<int>  triangleNeighbors;
  std::vector<unsigned> triangleNormalPoint;
  std::vector<double>  points;
  std::vector<float>  triangleSurfaces;
  std::vector<double> triangleCenterX;
  std::vector<double> triangleCenterY;
  std::vector<double>  betaCells;
  std::vector<unsigned>  claddingCellTypes;
  std::vector<float>  refractiveIndices;
  std::vector<float>  reflectivities;

  unsigned numberOfPoints = 0;
  unsigned numberOfTriangles = 0;
  unsigned numberOfLevels = 0;
  unsigned claddingNumber;
  float thickness = 1;
  float nTot = 0;
  float crystalTFluo = 0;
  double claddingAbsorption = 0;

  // Parse experimentdata from files
  if(fileToVector(rootPath + "triangleNormalPoint.txt", &triangleNormalPoint)) return meshs;
  if(fileToVector(rootPath + "betaVolume.txt", &betaVolume)) return meshs;
  if(fileToVector(rootPath + "forbiddenEdge.txt", &forbiddenEdge)) return meshs;
  if(fileToVector(rootPath + "triangleNeighbors.txt", &triangleNeighbors)) return meshs;
  if(fileToVector(rootPath + "triangleNormalsX.txt", &triangleNormalsX)) return meshs;
  if(fileToVector(rootPath + "triangleNormalsY.txt", &triangleNormalsY)) return meshs;
  if(fileToVector(rootPath + "triangleCenterX.txt", &triangleCenterX)) return meshs;
  if(fileToVector(rootPath + "triangleCenterY.txt", &triangleCenterY)) return meshs;
  if(fileToVector(rootPath + "points.txt", &points)) return meshs;
  if(fileToVector(rootPath + "trianglePointIndices.txt", &trianglePointIndices)) return meshs;
  if(fileToVector(rootPath + "triangleSurfaces.txt", &triangleSurfaces)) return meshs;
  if(fileToValue(rootPath  + "numberOfPoints.txt", numberOfPoints)) return meshs;
  if(fileToValue(rootPath  + "numberOfTriangles.txt", numberOfTriangles)) return meshs;
  if(fileToValue(rootPath  + "numberOfLevels.txt", numberOfLevels)) return meshs;
  if(fileToValue(rootPath  + "thickness.txt", thickness)) return meshs;
  if(fileToValue(rootPath  + "nTot.txt", nTot)) return meshs;
  if(fileToValue(rootPath  + "crystalTFluo.txt", crystalTFluo)) return meshs;
  if(fileToValue(rootPath  + "claddingNumber.txt", claddingNumber)) return meshs;
  if(fileToValue(rootPath  + "claddingAbsorption.txt", claddingAbsorption)) return meshs;
  if(fileToVector(rootPath + "betaCells.txt", &betaCells)) return meshs;
  if(fileToVector(rootPath + "claddingCellTypes.txt", &claddingCellTypes)) return meshs;
  if(fileToVector(rootPath + "refractiveIndices.txt", &refractiveIndices)) return meshs;
  if(fileToVector(rootPath + "reflectivities.txt", &reflectivities)) return meshs;

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

  for( unsigned i=0; i < maxGpus; i++){
    CUDA_CHECK_RETURN(cudaSetDevice(devices.at(i)) );
    meshs.push_back(
		    createMesh(trianglePointIndices,
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
			       claddingAbsorption
			       )
		    );
    cudaDeviceSynchronize();
  }

  return meshs;

}
