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

#include <any>
#include <filesystem>
#include <optional>
#include <ranges>
namespace fs = std::filesystem;


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
    )
    {
        std::stringstream ss;
        ss << CompSwitch::parallel_mode << " must be one of \""
           << ParallelMode::THREADED << "\"";
#if !defined(DISABLE_MPI)
        ss << ", \"" << ParallelMode::MPI << "\"";
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


void checkPositive(int i, const std::string& name, int threshhold=0){
    if(i < threshhold){
        verbosity |= V_ERROR;
        dout(V_ERROR) << name << " must have a positive argument!" << std::endl;
        throw std::invalid_argument(std::to_string(i));
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
    unsigned numberOfSamples=host_mesh.numberOfPoints * host_mesh.numberOfLevels;
    validateSampleRange(compute,numberOfSamples);
    auto wave=calculateSigmas(experiment.sigmaA,experiment.sigmaE,experiment.lambdaA,experiment.lambdaE,experiment.spectral);
    experiment.sigmaA=wave.sigmaAInterpolated;
    experiment.sigmaE=wave.sigmaEInterpolated;
    experiment.maxSigmaA=wave.maxSigmaA;
    experiment.maxSigmaE=wave.maxSigmaE;
    compute.devices = devices;
    compute.gpu_i=devices.front();

    mesh=std::move(createMeshes(host_mesh,devices));

    initializeResult(result, numberOfSamples);
    return 0;
}

int parse( const int argc,
        char** argv,
        ExperimentParameters& experiment,
        ComputeParameters& compute,
        std::vector<Mesh>& meshs,
        Result& result) {

    auto vm = parseCommandLine(argc, argv);
    printCommandLine(vm);

    // Set/Test device to run experiment with
    //TODO: this call takes a LOT of time (2-5s). Can this be avoided?
    //      maybe move this to a place where GPUs are actually needed (for_loops_clad doesn't even need GPUs!)
    auto &nr_gpuOpt=vm[CompSwitch::ngpus];
    checkPositive(nr_gpuOpt.as<int>(),nr_gpuOpt.name);
    std::vector<unsigned>devices = getFreeDevices(nr_gpuOpt.as<int>());
    ensureOrThrow(!devices.empty()," No cuda capable device found!");
    HostMesh host_mesh = createHostMeshFromFile(vm[ExpSwitch::input_path].as<fs::path>());
    validateHostMesh(host_mesh);
    meshs=std::move(createMeshes(host_mesh,devices));

    WavelengthData waveD = calculateSigmas(
            vm[ExpSwitch::input_path].as<fs::path>(),
            vm.isExplicitlySet(ExpSwitch::spectral) ? vm[ExpSwitch::spectral].as<int>() : 0
            );
    auto &min_rayOpt=vm[ExpSwitch::min_rays];
    auto &max_rayOpt=vm[ExpSwitch::max_rays];
    checkPositive(min_rayOpt.as<int>(),min_rayOpt.name);
    checkPositive(max_rayOpt.as<int>(),max_rayOpt.name);
    experiment = ExperimentParameters (
            static_cast<unsigned>(min_rayOpt.as<int>()),
            static_cast<unsigned>(max_rayOpt.as<int>()),
            waveD.sigmaAInterpolated,
            waveD.sigmaEInterpolated,
            waveD.maxSigmaA,
            waveD.maxSigmaE,
            vm[ExpSwitch::mse].as<double>(),
            vm[ExpSwitch::reflection].as<bool>() );

    auto &repetionOpt=vm[ExpSwitch::min_rays];
    auto &adaptiveStepOpt=vm[ExpSwitch::max_rays];
    checkPositive(repetionOpt.as<int>(),repetionOpt.name);
    checkPositive(adaptiveStepOpt.as<int>(),adaptiveStepOpt.name);
    compute = ComputeParameters (
    static_cast<unsigned>(repetionOpt.as<int>()),
    static_cast<unsigned>(adaptiveStepOpt.as<int>()),
            devices.at(0),
            vm[CompSwitch::device_mode].as<std::string>(),
            vm[CompSwitch::parallel_mode].as<std::string>(),
            vm[CompSwitch::write_vtk].as<bool>(),
            vm[ExpSwitch::input_path].as<fs::path>(),
            vm[ExpSwitch::output_path].as<fs::path>(),
            devices,
            vm[CompSwitch::min_sample_i].as<unsigned>(),
            vm[CompSwitch::max_sample_i].as<unsigned>(),
            nr_gpuOpt.as<int>());
    validateInput(experiment,compute,devices.size());
    validateSampleRange(compute,meshs[0].numberOfSamples);
    initializeResult(result, meshs[0].numberOfSamples);
    return 0;
}

inline CmdOptionsMap parseCommandLine(int argc, char** argv)
{
    CmdOptionsMap cfg{};

    for (int i = 1; i < argc; ++i)
    {
        std::string_view arg = argv[i];

        // allow --name=value form
        auto splitEq = [&](std::string_view token) -> std::pair<std::string_view, std::optional<std::string_view>> {
            auto pos = token.find('=');
            if (pos == std::string_view::npos) {
                return {token, std::nullopt};
            }
            return {token.substr(0, pos), token.substr(pos + 1)};
        };

        auto [name, inlineValue] = splitEq(arg);

        auto getValue = [&](std::string_view optName) -> std::string_view {
            if (inlineValue.has_value()) {
                return *inlineValue;
            }
            return requireValue(i, argc, argv, optName);
        };
        auto removePreDashes = [&](std::string_view optName) -> std::string_view
        {
            if (optName.empty()) {
                throw std::runtime_error("Invalid option name: empty string");
            }

            if (optName[0] != '-') {
                throw std::runtime_error(
                    "Invalid option name '" + std::string(optName) + "': expected leading '-' or '--'"
                );
            }

            const auto pos = optName.find_first_not_of('-');
            if (pos == std::string_view::npos) {
                throw std::runtime_error(
                    "Invalid option name '" + std::string(optName) + "': missing option name after dashes"
                );
            }

            return optName.substr(pos);
        };
        auto cleanName=removePreDashes(name);
        if (cleanName == "help" || cleanName == "h")
        {
            cfg.printDescription();
            std::terminate();
        }
        cfg.set(std::string(cleanName),std::string(getValue(name)));

    }
    cfg.checkRequired();
    return cfg;
}
void printCommandLine(const CmdOptionsMap& vm)
{
    for (const auto& [name, opt] : vm.options) {
        const std::string& shownValue = opt.hasValue ? opt.value : opt.defaultValue;
        dout(V_INFO) << name << ": " << shownValue << std::endl;
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

   return { cladAbsorption,
         static_cast<float>(totalSurface), //this cast might needlessly cause precision loss
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
         positionsOfNormalVectors};

}
