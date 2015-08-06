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
#include <sstream> /* stringstream */
#include <functional> /* bind, placeholders */
#define _USE_MATH_DEFINES /* make M_PI known for windows */
#include <cmath> /* M_PI */

#include <logging.hpp> 
#include <mesh.hpp>
#include <parser.hpp>
#include <interpolation.hpp> /* interpolateWavelength*/
#include <cuda_utils.hpp> /* getFreeDevices */


// includes for commandline parsing
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/variables_map.hpp>

#include <boost/filesystem.hpp> /* fs::path */
namespace fs = boost::filesystem;
namespace po = boost::program_options;


struct WavelengthData{
    std::vector<double>sigmaAInterpolated;
    std::vector<double>sigmaEInterpolated;
    double maxSigmaA = 0;
    double maxSigmaE = 0;
};


WavelengthData calculateSigmas(
        fs::path inputPath,
        unsigned lambdaResolution
        ){
        
    // Parse wavelengths from files
    std::vector<double> sigmaA  = fileToVector<double>(inputPath / "sigmaA.txt");
    std::vector<double> sigmaE  = fileToVector<double>(inputPath / "sigmaE.txt");
    std::vector<double> lambdaA = fileToVector<double>(inputPath / "lambdaA.txt");
    std::vector<double> lambdaE = fileToVector<double>(inputPath / "lambdaE.txt");

    lambdaResolution = std::max(lambdaResolution, (unsigned) lambdaA.size());
    lambdaResolution = std::max(lambdaResolution, (unsigned) lambdaE.size());
  
    assert(sigmaA.size() == lambdaA.size());
    assert(sigmaE.size() == lambdaE.size());

    // Interpolate sigmaA / sigmaE function
    WavelengthData waveD;
    waveD.sigmaAInterpolated = interpolateLinear(sigmaA, lambdaA, lambdaResolution);
    waveD.sigmaEInterpolated = interpolateLinear(sigmaE, lambdaE, lambdaResolution);
    assert(waveD.sigmaAInterpolated.size() == waveD.sigmaEInterpolated.size());

    // Calc max sigmaA / sigmaE
    for(unsigned i = 0; i < sigmaE.size(); ++i){
        if(sigmaE.at(i) > waveD.maxSigmaE){
            waveD.maxSigmaE = sigmaE.at(i);
            waveD.maxSigmaA = sigmaA.at(i);
        }
    }
    return waveD;
}


void checkPositive(int i, std::string name){
    if(i < 0){
        verbosity |= V_ERROR;
        dout(V_ERROR) << name << " must have a positive argument!" << std::endl;
        throw po::invalid_option_value(std::to_string(i));
    }
}


int parse( const int argc,
        char** argv,
        ExperimentParameters& experiment,
        ComputeParameters& compute,
        std::vector<Mesh>& meshs,
        Result& result) {

    bool writeVtk = false;

    Modifiable_variables_map vm = parseCommandLine(argc, argv);

    printCommandLine(vm);

    // Set/Test device to run experiment with
    //TODO: this call takes a LOT of time (2-5s). Can this be avoided?
    //      maybe move this to a place where GPUs are actually needed (for_loops_clad doesn't even need GPUs!)
    std::vector<unsigned>devices = getFreeDevices(vm[CompSwitch::ngpus].as<int>());

    vm = checkParameterValidity(vm, devices.size());

    meshs = parseMesh(vm[ExpSwitch::input_path].as<fs::path>(), devices);

    vm = checkSampleRange(vm, meshs[0].numberOfSamples);

    WavelengthData waveD = calculateSigmas(
            vm[ExpSwitch::input_path].as<fs::path>(),
            vm.count(ExpSwitch::spectral) ? vm[ExpSwitch::spectral].as<int>() : 0
            );

    experiment = ExperimentParameters (
            vm[ExpSwitch::min_rays].as<int>(),
            vm[ExpSwitch::max_rays].as<int>(),
            waveD.sigmaAInterpolated,
            waveD.sigmaEInterpolated,
            waveD.maxSigmaA,
            waveD.maxSigmaE,
            vm[ExpSwitch::mse].as<double>(),
            vm[ExpSwitch::reflection].as<bool>() );


    compute = ComputeParameters (
            vm[CompSwitch::repetitions].as<int>(),
            vm[CompSwitch::adaptive_steps].as<int>(),
            devices.at(0),
            vm[CompSwitch::device_mode].as<std::string>(),
            vm[CompSwitch::parallel_mode].as<std::string>(),
            writeVtk,
            vm[ExpSwitch::input_path].as<fs::path>(),
            vm[ExpSwitch::output_path].as<fs::path>(),
            devices,
            vm[CompSwitch::min_sample_i].as<int>(),
            vm[CompSwitch::max_sample_i].as<int>());


    std::vector<float>    phiAse(meshs[0].numberOfSamples, 0);
    std::vector<double>   mse(meshs[0].numberOfSamples, 100000);
    std::vector<unsigned> totalRays(meshs[0].numberOfSamples, 0);
    std::vector<double>   dndtAse(meshs[0].numberOfSamples, 0);

    result = Result( phiAse,
            mse,
            totalRays,
            dndtAse );

    return 0;
}


po::variables_map parseCommandLine(const int argc, char** argv) {

    po::options_description experiment_options( "Experiment Options" );
    experiment_options.add_options()
        ( std::string(ExpSwitch::input_path + ",i").c_str(),
          po::value<fs::path> ()->required(),
          "The path to a folder that contains the input files")
        ( std::string(ExpSwitch::output_path + ",o").c_str(),
          po::value<fs::path> ()->required(),
          "The path to a folder that contains the output files")
        ( ExpSwitch::min_rays.c_str(),
          po::value<int> ()
          ->default_value(100000)
          ->notifier(std::bind(checkPositive, std::placeholders::_1, ExpSwitch::min_rays)),
          "The minimal number of rays to use for each sample point")
        ( ExpSwitch::max_rays.c_str(),
          po::value<int> ()
          ->default_value(100000)
          ->notifier(std::bind(checkPositive, std::placeholders::_1, ExpSwitch::min_rays)),
          "The maximal number of rays to use for each sample point")
        ( std::string(ExpSwitch::mse + ",m").c_str(),
          po::value<double> ()->default_value(0.1,"0.1"),
          "The MSE threshold used for adaptive/repetitive sampling")
        ( ExpSwitch::reflection.c_str(),
          po::value<bool> ()->default_value(true),
          "use reflections or not")
        ( ExpSwitch::spectral.c_str(),
          po::value<int> (),
          "The number of samples used to interpolate spectral intensities")
        ;

    po::options_description compute_options( "Compute Options" );
    compute_options.add_options()
        ( CompSwitch::parallel_mode.c_str(),
          po::value<std::string> ()->default_value("threaded"),
          std::string("Set the preferred way of parellelization (mpi, graybat, threaded), only valid with --"
              + CompSwitch::device_mode + "=gpu").c_str())
        ( CompSwitch::device_mode.c_str(),
          po::value<std::string> ()->default_value("gpu"),
          "Set the device to run the calculation (cpu, gpu)")
        ( std::string(CompSwitch::ngpus + ",g").c_str(),
          po::value<int> ()
          ->default_value(1)
          ->notifier(std::bind(checkPositive, std::placeholders::_1, ExpSwitch::min_rays)),
          std::string("The maximum number of GPUs to b used on a single node. Should be left unchanged for --"
              + CompSwitch::parallel_mode + "=mpi").c_str())
        ( std::string(CompSwitch::repetitions + ",r").c_str(),
          po::value<int> ()
          ->default_value(4)
          ->notifier(std::bind(checkPositive, std::placeholders::_1, ExpSwitch::min_rays)),
          "The number of repetitions to try, before the number of rays is increased by adaptive sampling")
        ( std::string(CompSwitch::adaptive_steps + ",a").c_str(),
          po::value<int> ()
          ->default_value(5)
          ->notifier(std::bind(checkPositive, std::placeholders::_1, ExpSwitch::min_rays)),
          std::string("The number of adaptive sampling steps that are used to split the range between "
              + ExpSwitch::min_rays + " and " + ExpSwitch::max_rays).c_str())
        ( CompSwitch::min_sample_i.c_str(),
          po::value<int> ()
          ->notifier(std::bind(checkPositive, std::placeholders::_1, ExpSwitch::min_rays)),
          "The the minimal index of sample points to simulate")
        ( CompSwitch::max_sample_i.c_str(),
          po::value<int> ()
          ->notifier(std::bind(checkPositive, std::placeholders::_1, ExpSwitch::min_rays)),
          "The the maximal index of sample points to simulate")
        ;

    po::options_description generic_options( "Generic Options" );
    generic_options.add_options()
        ( "verbosity,v",
          po::value<int> ()
          ->default_value(V_ERROR | V_INFO | V_WARNING | V_PROGRESS | V_STAT)
          ->notifier(std::bind(checkPositive, std::placeholders::_1, ExpSwitch::min_rays)),
          "Set the verbosity levels:\n\
          \t 0 (quiet)\n\
          \t 1 (error)\n\
          \t 2 (warning)\n\
          \t 4 (info)\n\
          \t 8 (statistics)\n\
          \t 16 (debug)\n\
          \t 32 (progressbar)\n")
        ;

    po::options_description cmd_only_options;
    cmd_only_options.add_options()
        ( "config,c",
          po::value<fs::path> (),
          "location of an optional config file")
        ( "help,h",
          "print this help message and exit" )
        ;

    po::variables_map vm;

    po::options_description cmdline_options;
    cmdline_options.add(experiment_options).add(compute_options).add(generic_options).add(cmd_only_options);
    po::store(po::parse_command_line( argc, argv, cmdline_options ), vm);

    if(vm.count("help")){
        verbosity |= V_NOLABEL;
        dout(V_NOLABEL) << "Usage: " << argv[0] << " -i <input folder> -o <output folder> [options] " << std::endl;
        dout(V_NOLABEL) << cmdline_options << std::endl;
        exit(0);
    }

    verbosity = vm["verbosity"].as<int>();

    if(vm.count("config")){
        fs::path configPath(vm["config"].as<fs::path>());
        if(fs::exists(configPath)){
            po::options_description configfile_options;
            configfile_options.add(experiment_options).add(compute_options).add(generic_options);
            std::ifstream configPathStream(configPath.string());
            po::store(po::parse_config_file( configPathStream, configfile_options ), vm);
        }else{
            dout(V_ERROR) << "Configuration file could not be read. Possible problems: path "
                << configPath << "does not exist or insufficient permissions" << std::endl;
            exit(1);
        }
    }

    po::notify(vm);
    return vm;
}


void printCommandLine(const Modifiable_variables_map vm){

    for (const auto& it : vm) {
        std::stringstream ss;
        auto& value = it.second.value();
        if      (auto v = boost::any_cast<uint32_t>   (&value)) ss << *v;
        else if (auto v = boost::any_cast<std::string>(&value)) ss << *v;
        else if (auto v = boost::any_cast<int>        (&value)) ss << *v;
        else if (auto v = boost::any_cast<bool>       (&value)) ss << *v;
        else if (auto v = boost::any_cast<fs::path>   (&value)) ss << *v;
        else if (auto v = boost::any_cast<double>     (&value)) ss << *v;
        else{
            dout(V_WARNING) << it.first << ": cast error" << std::endl;
            continue;
        }
        dout(V_INFO) << it.first << ": " << ss.str() << std::endl;
    }
}


Modifiable_variables_map checkParameterValidity(Modifiable_variables_map vm, const unsigned deviceCount){
        
    double mseThreshold = vm[ExpSwitch::mse].as<double>();
    unsigned maxRaysPerSample = vm[ExpSwitch::max_rays].as<int>();
    unsigned maxgpus = vm[CompSwitch::ngpus].as<int>();

    const unsigned minRaysPerSample = vm[ExpSwitch::min_rays].as<int>();
    const fs::path inputPath = vm[ExpSwitch::input_path].as<fs::path>();
    const fs::path outputPath = vm[ExpSwitch::output_path].as<fs::path>();
    const std::string deviceMode = vm[CompSwitch::device_mode].as<std::string>();
    const std::string parallelMode = vm[CompSwitch::parallel_mode].as<std::string>();
    const unsigned maxRepetitions = vm[CompSwitch::repetitions].as<int>();
    const unsigned adaptiveSteps = vm[CompSwitch::adaptive_steps].as<int>();

    if(deviceMode != DeviceMode::CPU && deviceMode != DeviceMode::GPU){
        dout(V_ERROR) << CompSwitch::device_mode << " must be either \""
            << DeviceMode::CPU << "\" or \"" << DeviceMode::GPU << "\"" << std::endl;
        exit(1);
    }

    if (deviceMode == DeviceMode::CPU && parallelMode != ParallelMode::THREADED){
        dout(V_WARNING) << CompSwitch::device_mode << " \"" << DeviceMode::CPU << "\" does only support "
            << CompSwitch::parallel_mode << " \"threaded\"! (will be ignored)" << std::endl;
    }

    if (parallelMode != ParallelMode::GRAYBAT
            && parallelMode != ParallelMode::MPI
            && parallelMode != ParallelMode::THREADED) {
        dout(V_ERROR) << CompSwitch::parallel_mode << " must be either \"" << ParallelMode::MPI << "\", \""
            << ParallelMode::GRAYBAT << "\", or \"" << ParallelMode::THREADED << "\"" << std::endl;
        exit(1);
    }

    if (minRaysPerSample == 0) {
        dout(V_ERROR) << "Please specify the number of rays per sample Point with --" << ExpSwitch::min_rays << "=RAYS" << std::endl;
        exit(1);
    }

    if (!exists(inputPath) || !is_directory(inputPath)){
        dout(V_ERROR) << "The specified input path does not exist, is no directory, or has insufficient permissions." <<
            " Please specify a correct path by --" << ExpSwitch::output_path << "=[path]" << std::endl;
        exit(1);
    }else{
        if(is_empty(inputPath)) {
            dout(V_ERROR) << "The specified input folder " << ExpSwitch::input_path << " is empty!" << std::endl;
            exit(1);
        }
    }

    if (!exists(outputPath) || !is_directory(outputPath)) {
        dout(V_ERROR) << "The specified output path does not exist (or permission denied)."
            << " Please specify a correct folder by --" << ExpSwitch::output_path << "=[path]" << std::endl;
        exit(1);
    }

    if(maxRaysPerSample < minRaysPerSample){
        dout(V_WARNING) << ExpSwitch::max_rays << " < " << ExpSwitch::min_rays << ". Auto-increasing "
            << ExpSwitch::max_rays << " (will be non-adaptive!)" << std::endl;
        maxRaysPerSample = minRaysPerSample;
    }

    if(maxgpus > deviceCount){
        dout(V_ERROR) << "You don't have so many devices, use --" << CompSwitch::ngpus << "=" << deviceCount << std::endl;
        exit(1);
    }

    if(maxgpus == 0){
        maxgpus = deviceCount;
    }

    if( vm.count(CompSwitch::min_sample_i) + vm.count(CompSwitch::max_sample_i) == 2){
        const unsigned minSampleRange = vm[CompSwitch::min_sample_i].as<int>();
        const unsigned maxSampleRange = vm[CompSwitch::max_sample_i].as<int>();

        if(maxSampleRange < minSampleRange){
            dout(V_ERROR) << CompSwitch::max_sample_i << " < " << CompSwitch::min_sample_i << "!" << std::endl;
            exit(1);
        }

        unsigned samplesForNode = maxSampleRange-minSampleRange+1;
        if(maxgpus > samplesForNode){
            dout(V_WARNING) << "More GPUs requested than there are sample points. Number of used GPUs reduced to " << samplesForNode << std::endl;
            maxgpus = samplesForNode;
        }
    }

    if(verbosity >= 64){
        verbosity = 63;
        dout(V_WARNING) << "Verbosity level should be between 0 (quiet) and 63 (all). Levels can be bitmasked together." << std::endl;
    }
    if(maxRepetitions < 1){
        dout(V_ERROR) << "At least 1 repetition is necessary!" << std::endl;
    }

    if(adaptiveSteps < 1){
        dout(V_ERROR) << "At least 1 adaptive step is necessary!" << std::endl;
    }

    if(mseThreshold == 0){
        mseThreshold = 1000.;
    }

    vm[ExpSwitch::max_rays].value() = boost::any(static_cast<int>(maxRaysPerSample));
    vm[CompSwitch::ngpus].value() = boost::any(static_cast<int>(maxgpus));
    vm[ExpSwitch::mse].value() = boost::any(static_cast<double>(mseThreshold));

    return vm;
}


Modifiable_variables_map checkSampleRange(Modifiable_variables_map vm, const unsigned numberOfSamples){
    unsigned minCount = vm.count(CompSwitch::min_sample_i);
    unsigned maxCount = vm.count(CompSwitch::max_sample_i);

    if(minCount+maxCount < 1){
        dout(V_WARNING) << CompSwitch::min_sample_i << "/" << CompSwitch::max_sample_i
            << " not set! Assuming a sampling point range of " << "0 to " << numberOfSamples-1 << std::endl;
        vm[CompSwitch::min_sample_i].value() = boost::any(0);
        vm[CompSwitch::max_sample_i].value() = boost::any(numberOfSamples-1);
        return vm;
    }

    if(minCount+maxCount < 2){
        dout(V_ERROR) << "Options " << CompSwitch::min_sample_i << "/" << CompSwitch::max_sample_i
            << " must be used together! (Allowed Range from 0 to " << numberOfSamples-1 << ")";
        exit(1);
    }

    // if the code did not terminate, both min and max are defined inside the map
    if(static_cast<unsigned>(vm[CompSwitch::min_sample_i].as<int>()) > numberOfSamples){
        dout(V_ERROR) << CompSwitch::min_sample_i << " is out of range! (There are only " << numberOfSamples << " sampling points)";
        exit(1);
    }

    if(static_cast<unsigned>(vm[CompSwitch::max_sample_i].as<int>()) > numberOfSamples){
        dout(V_ERROR) << CompSwitch::max_sample_i << " is out of range! (There are only " << numberOfSamples << " sampling points)";
        exit(1);
    }

    return vm;
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
			    std::vector<unsigned> devices) {

  std::vector<Mesh> meshs;

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

  for( unsigned i=0; i < devices.size(); i++){
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
