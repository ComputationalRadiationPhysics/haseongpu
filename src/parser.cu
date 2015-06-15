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



#include <string> /* string */
#include <vector> /* vector */
#include <assert.h>
#include <stdlib.h> /* exit() */

#include <logging.hpp> 
#include <mesh.hpp>
#include <parser.hpp>
#include <interpolation.hpp> /* interpolateWavelength*/


// includes for commandline parsing
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/variables_map.hpp>

#include <boost/filesystem.hpp> /* fs::path */
namespace fs = boost::filesystem;


int parse( const int argc,
	   char** argv,
	   ExperimentParameters& experiment,
	   ComputeParameters& compute,
	   std::vector<Mesh>& meshs,
	   Result& result) {


    unsigned minRaysPerSample = 0;
    unsigned maxRaysPerSample = 0;
    unsigned maxRepetitions = 4;
    unsigned lambdaResolution = 0;
    std::string compareLocation("");
    bool writeVtk = false;
    bool useReflections = false;
    std::vector<unsigned> devices; 
    unsigned maxGpus = 0;
    DeviceMode deviceMode = NO_DEVICE_MODE;
    ParallelMode parallelMode = NO_PARALLEL_MODE;
    int minSampleRange = -1;
    int maxSampleRange = -1;
    unsigned usedGpus  = 0;
    fs::path inputPath;
    fs::path outputPath;
    double mseThreshold = 0;

    // Wavelength data
    std::vector<double> sigmaA;
    std::vector<double> sigmaE;
    std::vector<double> lambdaA;
    std::vector<double> lambdaE;

    // Parse Commandline

    // Can the following block be moved into
    // a parse function with a thin interface ?

    //
    // BEGIN PARSE BLOCK
    //
    parseCommandLine(argc,
		     argv,
		     &minRaysPerSample, // exp
		     &maxRaysPerSample, // exp
		     &inputPath, // opt -
		     &writeVtk,  // opt -
		     &deviceMode,// opt -
		     &parallelMode, // opt -
		     &useReflections, // exp
		     &maxGpus, // opt -
		     &minSampleRange, // opt -
		     &maxSampleRange, // opt -
		     &maxRepetitions, // exp
		     &outputPath, // opt -
		     &mseThreshold, // exp
		     &lambdaResolution); // opt

    printCommandLine(minRaysPerSample,
		     maxRaysPerSample,
		     inputPath,
		     writeVtk,
		     compareLocation,
		     deviceMode,
		     parallelMode,
		     useReflections,
		     maxGpus,
		     minSampleRange,
		     maxSampleRange,
		     maxRepetitions,
		     outputPath,
		     mseThreshold);
  
    // Set/Test device to run experiment with
    //
    //TODO: this call takes a LOT of time (2-5s). Can this be avoided?
    //TODO: maybe move this to a place where GPUs are actually needed (for_loops_clad doesn't even need GPUs!)
    devices = getFreeDevices(maxGpus);

    // sanity checks
    if(checkParameterValidity(argc,
			      minRaysPerSample,
			      &maxRaysPerSample,
			      inputPath,
			      devices.size(),
			      deviceMode,
			      parallelMode,
			      &maxGpus,
			      minSampleRange,
			      maxSampleRange,
			      maxRepetitions,
			      outputPath,
			      &mseThreshold)) return 1;

  
    dout(V_INFO) << "parameter validity was checked!" << std::endl;

    // Parse wavelengths from files
    if(fileToVector(inputPath, "sigmaA.txt",  &sigmaA))   return 1;
    if(fileToVector(inputPath, "sigmaE.txt",  &sigmaE))   return 1;
    if(fileToVector(inputPath, "lambdaA.txt", &lambdaA)) return 1;
    if(fileToVector(inputPath, "lambdaE.txt", &lambdaE)) return 1;
    lambdaResolution = std::max(lambdaResolution, (unsigned) lambdaA.size());
    lambdaResolution = std::max(lambdaResolution, (unsigned) lambdaE.size());
  
    assert(sigmaA.size() == lambdaA.size());
    assert(sigmaE.size() == lambdaE.size());

    // Interpolate sigmaA / sigmaE function
    std::vector<double> sigmaAInterpolated = interpolateLinear(sigmaA, lambdaA, lambdaResolution); // exp
    std::vector<double> sigmaEInterpolated = interpolateLinear(sigmaE, lambdaE, lambdaResolution); // exp
    assert(sigmaAInterpolated.size() == sigmaEInterpolated.size());

    // Calc max sigmaA / sigmaE
    // TODO replace by some algorithm
    double maxSigmaE = 0.0;
    double maxSigmaA = 0.0;
    for(unsigned i = 0; i < sigmaE.size(); ++i){
	if(sigmaE.at(i) > maxSigmaE){
	    maxSigmaE = sigmaE.at(i);
	    maxSigmaA = sigmaA.at(i);
	}
    }

    // Parse experientdata and fill mesh
    meshs = parseMesh(inputPath, devices, maxGpus);

    checkSampleRange(&minSampleRange,&maxSampleRange,meshs[0].numberOfSamples);

    // Solution vector
    std::vector<double>   dndtAse(meshs[0].numberOfSamples, 0); 
    std::vector<float>    phiAse(meshs[0].numberOfSamples, 0); // exp
    std::vector<double>   mse(meshs[0].numberOfSamples, 100000); // exp
    std::vector<unsigned> totalRays(meshs[0].numberOfSamples, 0); // exp

    //
    // END PARSE BLOCK
    //

    experiment = ExperimentParameters ( minRaysPerSample,
					maxRaysPerSample,
					sigmaAInterpolated,
					sigmaEInterpolated,
					maxSigmaA,
					maxSigmaE,
					mseThreshold,
					useReflections );

    compute = ComputeParameters ( maxRepetitions,
				  devices.at(0),
				  deviceMode,
				  parallelMode,
				  writeVtk,
				  inputPath,
				  outputPath,
				  devices,
				  minSampleRange,
				  maxSampleRange);

    result = Result( phiAse, 
		     mse, 
		     totalRays,
		     dndtAse );

    return 0;

}
    



void parseCommandLine(
    const int argc,
    char** argv,
    unsigned *raysPerSample,
    unsigned *maxRaysPerSample,
    fs::path *inputPath,
    bool *writeVtk,
    DeviceMode *deviceMode,
    ParallelMode *parallelMode,
    bool *useReflections,
    unsigned *maxgpus,
    int *minSample_i,
    int *maxSample_i,
    unsigned *maxRepetitions,
    fs::path *outputPath,
    double *mseThreshold,
    unsigned *lambdaResolution
    ) {

  std::string dMode;
  std::string pMode;
  namespace po = boost::program_options;
  po::options_description desc( "Allowed options" );
  desc.add_options()
    ( "help,h",
      "print this help message and exit" )
    ( "verbosity,v",
      po::value<unsigned> ( &verbosity )->default_value(V_ERROR | V_INFO | V_WARNING | V_PROGRESS | V_STAT),
      "Set the verbosity levels:\n\
      \t 0 (quiet)\n\
      \t 1 (error)\n\
      \t 2 (warning)\n\
      \t 4 (info)\n\
      \t 8 (statistics)\n\
      \t 16 (debug)\n\
      \t 32 (progressbar)\n")
    ( "device-mode", po::value<std::string> (&dMode)->default_value("gpu"),
      "Set the device to run the calculation (cpu, gpu)")
    ( "parallel-mode", po::value<std::string> (&pMode)->default_value("threaded"),
      "Set the preferred way of parellelization (mpi, threaded), only valid with device-mode=gpu")
    ( "reflection", po::value<bool> (useReflections)->default_value(true),
      "use reflections or not")
    ( "min-rays", po::value<unsigned> (raysPerSample)->default_value(100000),
      "The minimal number of rays to use for each sample point")
    ( "max-rays", po::value<unsigned> (maxRaysPerSample)->default_value(100000),
      "The maximal number of rays to use for each sample point")
    ( "input-path,i", po::value<fs::path> (inputPath)->default_value(fs::path("input/")),
      "The path to a folder that contains the input files")
    ( "output-path,o", po::value<fs::path> (outputPath)->default_value(fs::path("output/")),
      "The path to a folder that contains the output files")
    ( "ngpus,g", po::value<unsigned> (maxgpus)->default_value(1),
      "The maximum number of GPUs to b used on a single node. Should be left unchanged for --parallel-mode=mpi")
    ( "min-sample-i", po::value<int> (minSample_i),
      "The the minimal index of sample points to simulate")
    ( "max-sample-i", po::value<int> (maxSample_i),
      "The the maximal index of sample points to simulate")
    ( "mse-threshold,m", po::value<double> (mseThreshold)->default_value(0.1,"0.1"),
      "The MSE threshold used for adaptive/repetitive sampling")
    ( "spectral-resolution", po::value<unsigned> (lambdaResolution),
      "The number of samples used to interpolate spectral intensities")
    ( "repetitions,r", po::value<unsigned> (maxRepetitions)->default_value(4),
      "The number of repetitions to try, before the number of rays is increased by adaptive sampling");

  po::variables_map vm;
  po::store(po::parse_command_line( argc, argv, desc ), vm);
  po::notify(vm);

  if(vm.count("help")){
    std::cout << "Usage: " << argv[0] << " [options] " << std::endl;
    std::cout << std::endl;
    std::cout << desc << std::endl;
    exit(0);
  }

  if (pMode == "threaded")
    *parallelMode = THREADED_PARALLEL_MODE;
  else if (pMode == "mpi")
    *parallelMode = MPI_PARALLEL_MODE;
  else if (pMode == "graybat")
    *parallelMode = GRAYBAT_PARALLEL_MODE;
  else
    *parallelMode = NO_PARALLEL_MODE;

  if(dMode == "gpu")
    *deviceMode = GPU_DEVICE_MODE;
  else if (dMode == "cpu")
    *deviceMode = CPU_DEVICE_MODE;
  else
    *deviceMode = NO_DEVICE_MODE;
      
  // append trailing folder separator, if necessary
  *outputPath /= "";
  *inputPath /= "";



//    dout(V_ERROR) << "  --compare=<location of vtk-file to compare with>" << std::endl;
//    write-vtk

}
void printCommandLine(
    unsigned raysPerSample,
    unsigned maxRaysPerSample,
    const fs::path inputPath,
    bool writeVtk,
    const fs::path compareLocation,
    const DeviceMode dMode,
    const ParallelMode pMode,
    bool useReflections,
    unsigned maxgpus,
    int minSample_i,
    int maxSample_i,
    unsigned maxRepetitions,
    const fs::path outputPath,
    double mseThreshold){
    
  dout(V_INFO) << "raysPerSample: " << raysPerSample << std::endl;
  dout(V_INFO) << "maxRaysPerSample: " << maxRaysPerSample << std::endl;
  dout(V_INFO) << "inputPath: " << inputPath << std::endl;
  dout(V_INFO) << "outputPath: " << outputPath << std::endl;
  dout(V_INFO) << "device-mode: " << dMode << std::endl;
  dout(V_INFO) << "parallel-mode: " << pMode << std::endl;
  dout(V_INFO) << "useReflections: " << useReflections << std::endl;
  dout(V_INFO) << "maxgpus: " << maxgpus << std::endl;
  dout(V_INFO) << "minSample_i: " << minSample_i << std::endl;
  dout(V_INFO) << "maxSample_i:" << maxSample_i << std::endl;
  dout(V_INFO) << "maxRepetitions: " << maxRepetitions << std::endl;
  dout(V_INFO) << "mseThreshold: " << mseThreshold << std::endl;
}

int checkParameterValidity(
    const int argc,
    const unsigned raysPerSample,
    unsigned *maxRaysPerSample,
    const fs::path inputPath,
    const unsigned deviceCount,
    const DeviceMode deviceMode,
    const ParallelMode parallelMode,
    unsigned *maxgpus,
    const int minSampleRange,
    const int maxSampleRange,
    const unsigned maxRepetitions,
    const fs::path outputPath,
    double *mseThreshold
    ) {

  if (deviceMode == NO_DEVICE_MODE) {
    dout(V_ERROR) << "device-mode must be either \"gpu\" or \"cpu\" " << std::endl;
    return 1;
  }

  if (deviceMode == CPU_DEVICE_MODE && parallelMode == MPI_PARALLEL_MODE){
    dout(V_WARNING) << "device-mode \"cpu\" does not support parallel-mode \"mpi\"! (will be ignored)" << std::endl;
  }

  if (parallelMode == NO_PARALLEL_MODE) {
    dout(V_ERROR) << "parallel-mode must be either \"mpi\" or \"threaded\" " << std::endl;
    return 1;
  }

  if (raysPerSample == 0) {
    dout(V_ERROR) << "Please specify the number of rays per sample Point with --rays=RAYS" << std::endl;
    return 1;
  }

  if (!exists(inputPath) || !is_directory(inputPath) || is_empty(inputPath)) {
    dout(V_ERROR) << "Please specify the experiment's location with --input-path=PATH_TO_EXPERIMENT" << std::endl;
    return 1;
  }

  if (!exists(outputPath) || !is_directory(outputPath)) {
    dout(V_ERROR) << "Please specify the output location with --output-path=PATH_TO_EXPERIMENT" << std::endl;
    return 1;
  }


  if(*maxRaysPerSample < raysPerSample){
    dout(V_WARNING) << "maxRays < raysPerSample. Increasing maxRays to " << raysPerSample << " (will be non-adaptive!)" << std::endl;
    *maxRaysPerSample = raysPerSample;
  }

  if(*maxgpus > deviceCount){ dout(V_ERROR) << "You don't have so many devices, use --ngpus=" << deviceCount << std::endl;
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

  if((*minSampleRange < -1 || *minSampleRange >= (int)numberOfSamples)){
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
std::vector<Mesh> parseMesh(const fs::path rootPath,
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
  //fs::path triangleNormalPointPath(rootPath);
  if(fileToVector(rootPath, "triangleNormalPoint.txt", &triangleNormalPoint)) return meshs;
  if(fileToVector(rootPath, "betaVolume.txt", &betaVolume)) return meshs;
  if(fileToVector(rootPath, "forbiddenEdge.txt", &forbiddenEdge)) return meshs;
  if(fileToVector(rootPath, "triangleNeighbors.txt", &triangleNeighbors)) return meshs;
  if(fileToVector(rootPath, "triangleNormalsX.txt", &triangleNormalsX)) return meshs;
  if(fileToVector(rootPath, "triangleNormalsY.txt", &triangleNormalsY)) return meshs;
  if(fileToVector(rootPath, "triangleCenterX.txt", &triangleCenterX)) return meshs;
  if(fileToVector(rootPath, "triangleCenterY.txt", &triangleCenterY)) return meshs;
  if(fileToVector(rootPath, "points.txt", &points)) return meshs;
  if(fileToVector(rootPath, "trianglePointIndices.txt", &trianglePointIndices)) return meshs;
  if(fileToVector(rootPath, "triangleSurfaces.txt", &triangleSurfaces)) return meshs;
  if(fileToValue(rootPath, "numberOfPoints.txt", numberOfPoints)) return meshs;
  if(fileToValue(rootPath, "numberOfTriangles.txt", numberOfTriangles)) return meshs;
  if(fileToValue(rootPath, "numberOfLevels.txt", numberOfLevels)) return meshs;
  if(fileToValue(rootPath, "thickness.txt", thickness)) return meshs;
  if(fileToValue(rootPath, "nTot.txt", nTot)) return meshs;
  if(fileToValue(rootPath, "crystalTFluo.txt", crystalTFluo)) return meshs;
  if(fileToValue(rootPath, "claddingNumber.txt", claddingNumber)) return meshs;
  if(fileToValue(rootPath, "claddingAbsorption.txt", claddingAbsorption)) return meshs;
  if(fileToVector(rootPath, "betaCells.txt", &betaCells)) return meshs;
  if(fileToVector(rootPath, "claddingCellTypes.txt", &claddingCellTypes)) return meshs;
  if(fileToVector(rootPath, "refractiveIndices.txt", &refractiveIndices)) return meshs;
  if(fileToVector(rootPath, "reflectivities.txt", &reflectivities)) return meshs;

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
