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

#include <cmath> /* M_PI */
#include <string> /* string */
#include <vector> /* vector */
#include <cassert>
#include <algorithm>
#include <cstdlib> /* exit() */
#include <sstream> /* stringstream */
#include <functional> /* bind, placeholders */

#include <cuda_runtime_api.h> /* cudaSetDevice */

#include <cudachecks.hpp> /* CUDA_CHECK_RETURN */
#include <cuda_utils.hpp> /* getFreeDevices */
#include <types.hpp>
#include <logging.hpp>
#include <mesh.hpp>
#include <parser.hpp>
#include <interpolation.hpp> /* interpolateWavelength*/

// includes for commandline parsing
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/any.hpp> /* boost::any_cast */

#include <boost/filesystem/path.hpp> /* fs::path */
#include <boost/filesystem/operations.hpp>

namespace fs = boost::filesystem;
namespace po = boost::program_options;


struct WavelengthData{
    std::vector<double>sigmaAInterpolated;
    std::vector<double>sigmaEInterpolated;
    double maxSigmaA = 0;
    double maxSigmaE = 0;
};
/**
 *
 */
HostMesh createHostMeshFromFile(const fs::path rootPath)
{
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
    return HostMesh(
      trianglePointIndices,
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
  );
}

template <class T, class B, class E>
void assertRange([[maybe_unused]]const std::vector<T> &v,[[maybe_unused]] const B minElement,[[maybe_unused]]const E maxElement,[[maybe_unused]] const bool equals){
    if(equals){
        assert(*std::min_element(v.begin(),v.end()) == minElement);
        assert(*std::max_element(v.begin(),v.end()) == maxElement);
    }else{
        assert(*std::min_element(v.begin(),v.end()) >= minElement);
        assert(*std::max_element(v.begin(),v.end()) <= maxElement);
    }
}


template <class T, class B>
void assertMin([[maybe_unused]]const std::vector<T>& v,[[maybe_unused]]const  B minElement,const bool equals){
    if(equals){
        assert(*std::min_element(v.begin(),v.end()) == minElement);
    }else{
        assert(*std::min_element(v.begin(),v.end()) >= minElement);
    }
}

bool validateHostMesh(const HostMesh& mesh)
{
    // assert input sizes
    assert(mesh.numberOfPoints == (mesh.pointsVector.size() / 2));
    assert(mesh.numberOfTriangles == mesh.triangleIndices.size() / 3);
    assert(mesh.positionsOfNormalVectors.size() == mesh.numberOfTriangles * 3);
    assert(mesh.yOfTriangleCenter.size() == mesh.numberOfTriangles);
    assert(mesh.xOfTriangleCenter.size() == mesh.numberOfTriangles);
    assert(mesh.surfacesVector.size() == mesh.numberOfTriangles);
    assert(mesh.betaValuesVector.size() == mesh.numberOfTriangles * (mesh.numberOfLevels - 1));
    assert(mesh.xOfNormals.size() == mesh.numberOfTriangles * 3);
    assert(mesh.yOfNormals.size() == mesh.numberOfTriangles * 3);
    assert(mesh.triangleIndices.size() == mesh.numberOfTriangles * 3);
    assert(mesh.forbiddenVector.size() == mesh.numberOfTriangles * 3);
    assert(mesh.neighborsVector.size() == mesh.numberOfTriangles * 3);
    assert(mesh.betaCells.size() == mesh.numberOfPoints * mesh.numberOfLevels);
    assert(mesh.cellTypes.size() == mesh.numberOfTriangles);
    assert(mesh.refractiveIndices.size() == 4);
    assert(mesh.reflectivities.size() == (mesh.refractiveIndices.size() / 2) * mesh.numberOfTriangles);
    assert(mesh.cellTypes.size() == mesh.numberOfTriangles);

    // assert input data validity
    assertRange(mesh.positionsOfNormalVectors, 0u, unsigned(mesh.numberOfPoints - 1), true);
    assertMin(mesh.betaValuesVector, 0, false);
    assertRange(mesh.forbiddenVector, -1, 2, true);
    assertRange(mesh.neighborsVector, -1, int(mesh.numberOfTriangles - 1), true);
    assertRange(mesh.xOfNormals, -1, 1, false);
    assertRange(mesh.yOfNormals, -1, 1, false);

    assertRange(
        mesh.xOfTriangleCenter,
        *std::min_element(mesh.pointsVector.begin(), mesh.pointsVector.end()),
        *std::max_element(mesh.pointsVector.begin(), mesh.pointsVector.end()),
        false);

    assertRange(
        mesh.yOfTriangleCenter,
        *std::min_element(mesh.pointsVector.begin(), mesh.pointsVector.end()),
        *std::max_element(mesh.pointsVector.begin(), mesh.pointsVector.end()),
        false);

    assertRange(mesh.triangleIndices, 0u, unsigned(mesh.numberOfPoints - 1), true);
    assertMin(mesh.surfacesVector, 0, false);
    assertMin(mesh.betaCells, 0, false);
    assertRange(mesh.refractiveIndices, 0, 5, false);
    assertRange(mesh.reflectivities, 0, 1, false);

    return true;
}
void validateSampleRange(
    ComputeParameters& compute,
    const unsigned numberOfSamples)
{
    const bool minSet = compute.minSampleRange != std::numeric_limits<unsigned>::max();
    const bool maxSet = compute.maxSampleRange != std::numeric_limits<unsigned>::max();

    if(!minSet && !maxSet)
    {
        dout(V_WARNING) << CompSwitch::min_sample_i << "/" << CompSwitch::max_sample_i
                        << " not set! Assuming a sampling point range of 0 to "
                        << (numberOfSamples - 1) << std::endl;
        compute.minSampleRange = 0;
        compute.maxSampleRange = numberOfSamples - 1;
        return;
    }

    if(minSet != maxSet)
    {
        throw std::runtime_error(
            "Options " + std::string(CompSwitch::min_sample_i) + "/" +
            std::string(CompSwitch::max_sample_i) +
            " must be used together!");
    }

    if(compute.maxSampleRange < compute.minSampleRange)
    {
        throw std::runtime_error(
            std::string(CompSwitch::max_sample_i) + " < " +
            std::string(CompSwitch::min_sample_i) + "!");
    }

    if(compute.minSampleRange >= numberOfSamples)
    {
        throw std::runtime_error(
            std::string(CompSwitch::min_sample_i) +
            " is out of range! (There are only " +
            std::to_string(numberOfSamples) + " sampling points)");
    }

    if(compute.maxSampleRange >= numberOfSamples)
    {
        throw std::runtime_error(
            std::string(CompSwitch::max_sample_i) +
            " is out of range! (There are only " +
            std::to_string(numberOfSamples) + " sampling points)");
    }

    const unsigned samplesForNode = compute.maxSampleRange - compute.minSampleRange + 1;
    if(compute.maxGpus > samplesForNode)
    {
        dout(V_WARNING) << "More GPUs requested than there are sample points. "
                        << "Number of used GPUs reduced to " << samplesForNode << std::endl;
        compute.maxGpus = samplesForNode;
    }
}
inline void ensureOrThrow(bool cond, const std::string& msg)
{
    if(!cond)
        throw std::runtime_error(msg);
}
void validateInput(
    ExperimentParameters& experiment,
    ComputeParameters& compute,
    const unsigned deviceCount)
{
    if(compute.deviceMode != DeviceMode::CPU && compute.deviceMode != DeviceMode::GPU)
    {
        throw std::runtime_error(
            CompSwitch::device_mode + " must be either \"" +
            std::string(DeviceMode::CPU) + "\" or \"" +
            std::string(DeviceMode::GPU) + "\"");
    }

    if(compute.deviceMode == DeviceMode::CPU &&
       compute.parallelMode != ParallelMode::THREADED)
    {
        dout(V_WARNING) << CompSwitch::device_mode << " \"" << DeviceMode::CPU
                        << "\" does only support " << CompSwitch::parallel_mode
                        << " \"threaded\"! (will be ignored)" << std::endl;
    }

    if(compute.parallelMode != ParallelMode::THREADED
#if !defined(DISABLE_MPI)
       && compute.parallelMode != ParallelMode::MPI
#endif
#if defined(USE_GRAYBAT)
       && compute.parallelMode != ParallelMode::GRAYBAT
#endif
    )
    {
        std::stringstream ss;
        ss << CompSwitch::parallel_mode << " must be one of \""
           << ParallelMode::THREADED << "\"";
#if !defined(DISABLE_MPI)
        ss << ", \"" << ParallelMode::MPI << "\"";
#endif
#if defined(USE_GRAYBAT)
        ss << ", \"" << ParallelMode::GRAYBAT << "\"";
#endif
        throw std::runtime_error(ss.str());
    }

    ensureOrThrow(experiment.minRaysPerSample > 0,
        "Please specify a positive value for --" + ExpSwitch::min_rays);

    if(experiment.maxRaysPerSample < experiment.minRaysPerSample)
    {
        dout(V_WARNING) << ExpSwitch::max_rays << " < " << ExpSwitch::min_rays
                        << ". Auto-increasing " << ExpSwitch::max_rays
                        << " (will be non-adaptive!)" << std::endl;
        experiment.maxRaysPerSample = experiment.minRaysPerSample;
    }

    if(compute.maxGpus > deviceCount)
    {
        throw std::runtime_error(
            "You don't have so many devices, use --" +
            std::string(CompSwitch::ngpus) + "=" + std::to_string(deviceCount)+ " "+ std::to_string(compute.maxGpus));
    }

    if(compute.maxGpus == 0)
    {
        compute.maxGpus = deviceCount;
    }

    if(verbosity >= 64)
    {
        verbosity = 63;
        dout(V_WARNING) << "Verbosity level should be between 0 (quiet) and 63 (all). "
                        << "Levels can be bitmasked together." << std::endl;
    }

    ensureOrThrow(compute.maxRepetitions >= 1, "At least 1 repetition is necessary!");
    ensureOrThrow(compute.adaptiveSteps >= 1, "At least 1 adaptive step is necessary!");

    if(experiment.mseThreshold == 0.0)
    {
        experiment.mseThreshold = 1000.0;
    }

    if(!compute.inputPath.empty())
    {
        ensureOrThrow(fs::exists(compute.inputPath) && fs::is_directory(compute.inputPath),
            "The specified input path does not exist, is no directory, or has insufficient permissions.");

        ensureOrThrow(!fs::is_empty(compute.inputPath),
            "The specified input folder is empty.");
    }

    if(!compute.outputPath.empty())
    {
        ensureOrThrow(fs::exists(compute.outputPath) && fs::is_directory(compute.outputPath),
            "The specified output path does not exist (or permission denied).");
    }
}
WavelengthData calculateSigmas(std::vector<double> &sigmaA,std::vector<double> &sigmaE,std::vector<double> &lambdaA,std::vector<double> &lambdaE,unsigned lambdaResolution)
{
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

WavelengthData calculateSigmas(
        fs::path inputPath,
        unsigned lambdaResolution
        ){

    // Parse wavelengths from files
    std::vector<double> sigmaA  = fileToVector<double>(inputPath / "sigmaA.txt");
    std::vector<double> sigmaE  = fileToVector<double>(inputPath / "sigmaE.txt");
    std::vector<double> lambdaA = fileToVector<double>(inputPath / "lambdaA.txt");
    std::vector<double> lambdaE = fileToVector<double>(inputPath / "lambdaE.txt");
    return calculateSigmas(sigmaA, sigmaE, lambdaA, lambdaE, lambdaResolution);

}


void checkPositive(int i, const std::string& name){
    if(i < 0){
        verbosity |= V_ERROR;
        dout(V_ERROR) << name << " must have a positive argument!" << std::endl;
        throw po::invalid_option_value(std::to_string(i));
    }
}


void initializeResult(Result& result, const unsigned numberOfSamples)
{
    result = Result(
        std::vector<float>(numberOfSamples, 0.0f),
        std::vector<double>(numberOfSamples, 100000.0),
        std::vector<unsigned>(numberOfSamples, 0u),
        std::vector<double>(numberOfSamples, 0.0));
}
std::vector<Mesh> createMeshes(HostMesh& host_mesh,std::vector<unsigned> &devices)
{
    std::vector<Mesh> meshes;
    for(auto i:devices)
    {
        CUDA_CHECK_RETURN(cudaSetDevice(i) );
        meshes.emplace_back(host_mesh.toMesh());
        cudaDeviceSynchronize();

    }
    return meshes;
}
int pythonParse(ExperimentParameters& experiment,
    ComputeParameters& compute,
    HostMesh & host_mesh,
    std::vector<Mesh>& mesh,
    Result& result)
{
    std::vector<unsigned>devices = getFreeDevices(compute.maxGpus);
    ensureOrThrow(!devices.empty()," No cuda capable device found!");
    validateHostMesh(host_mesh);
    validateInput(experiment, compute, devices.size());
    validateSampleRange(compute,host_mesh.numberOfPoints * host_mesh.numberOfLevels);
    auto wave=calculateSigmas(experiment.sigmaA,experiment.sigmaE,experiment.lambdaA,experiment.lambdaE,experiment.spectral);
    experiment.sigmaA=wave.sigmaAInterpolated;
    experiment.sigmaE=wave.sigmaEInterpolated;
    experiment.maxSigmaA=wave.maxSigmaA;
    experiment.maxSigmaE=wave.maxSigmaE;
    compute.devices = devices;
    compute.gpu_i=devices.front();

    mesh=std::move(createMeshes(host_mesh,devices));

    initializeResult(result, host_mesh.numberOfPoints);
    return 0;
}

int parse( const int argc,
        char** argv,
        ExperimentParameters& experiment,
        ComputeParameters& compute,
        std::vector<Mesh>& meshs,
        Result& result) {

    Modifiable_variables_map vm = parseCommandLine(argc, argv);
    printCommandLine(vm);

    // Set/Test device to run experiment with
    //TODO: this call takes a LOT of time (2-5s). Can this be avoided?
    //      maybe move this to a place where GPUs are actually needed (for_loops_clad doesn't even need GPUs!)
    std::vector<unsigned>devices = getFreeDevices(vm[CompSwitch::ngpus].as<int>());
    ensureOrThrow(!devices.empty()," No cuda capable device found!");
    HostMesh host_mesh = createHostMeshFromFile(vm[ExpSwitch::input_path].as<fs::path>());
    validateHostMesh(host_mesh);
    meshs=std::move(createMeshes(host_mesh,devices));

    WavelengthData waveD = calculateSigmas(
            vm[ExpSwitch::input_path].as<fs::path>(),
            vm.count(ExpSwitch::spectral) ? vm[ExpSwitch::spectral].as<int>() : 0
            );

    experiment = ExperimentParameters (
            static_cast<unsigned>(vm[ExpSwitch::min_rays].as<int>()),
            static_cast<unsigned>(vm[ExpSwitch::max_rays].as<int>()),
            waveD.sigmaAInterpolated,
            waveD.sigmaEInterpolated,
            waveD.maxSigmaA,
            waveD.maxSigmaE,
            vm[ExpSwitch::mse].as<double>(),
            vm[ExpSwitch::reflection].as<bool>() );


    compute = ComputeParameters (
    static_cast<unsigned>(vm[CompSwitch::repetitions].as<int>()),
    static_cast<unsigned>(vm[CompSwitch::adaptive_steps].as<int>()),
            devices.at(0),
            vm[CompSwitch::device_mode].as<std::string>(),
            vm[CompSwitch::parallel_mode].as<std::string>(),
            vm[CompSwitch::write_vtk].as<bool>(),
            vm[ExpSwitch::input_path].as<fs::path>(),
            vm[ExpSwitch::output_path].as<fs::path>(),
            devices,
            vm[CompSwitch::min_sample_i].as<unsigned>(),
            vm[CompSwitch::max_sample_i].as<unsigned>(),
            vm[CompSwitch::ngpus].as<int>());
    validateInput(experiment,compute,devices.size());
    validateSampleRange(compute,meshs[0].numberOfSamples);
    initializeResult(result, meshs[0].numberOfSamples);
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
          ->notifier(std::bind(checkPositive, std::placeholders::_1, ExpSwitch::max_rays)),
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
          ->notifier(std::bind(checkPositive, std::placeholders::_1, CompSwitch::ngpus)),
          std::string("The maximum number of GPUs to b used on a single node. Should be left unchanged for --"
              + CompSwitch::parallel_mode + "=mpi").c_str())
        ( std::string(CompSwitch::repetitions + ",r").c_str(),
          po::value<int> ()
          ->default_value(4)
          ->notifier(std::bind(checkPositive, std::placeholders::_1, CompSwitch::repetitions)),
          "The number of repetitions to try, before the number of rays is increased by adaptive sampling")
        ( std::string(CompSwitch::adaptive_steps + ",a").c_str(),
          po::value<int> ()
          ->default_value(5)
          ->notifier(std::bind(checkPositive, std::placeholders::_1, CompSwitch::adaptive_steps)),
          std::string("The number of adaptive sampling steps that are used to split the range between "
              + ExpSwitch::min_rays + " and " + ExpSwitch::max_rays).c_str())
        ( CompSwitch::min_sample_i.c_str(),
          po::value<unsigned> ()
          ->default_value(-1),
          "The the minimum index of sample points to simulate")
        ( CompSwitch::max_sample_i.c_str(),
          po::value<unsigned> ()
          ->default_value(-1),
          "The the maximal index of sample points to simulate")
        ( CompSwitch::write_vtk.c_str(),
          po::value<bool> ()
          ->default_value(false),
          "Write VTK files of the computed ASE values")
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
        dout(V_INFO) << "[DEBUG] option '" << it.first
                     << "' stored type: " << value.type().name() << std::endl;

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
