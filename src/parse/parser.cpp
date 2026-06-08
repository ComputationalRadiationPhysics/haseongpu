/**
 * Copyright 2015 Erik Zenker, Carlchristian Eckert, Marius Melzer
 * Copyright 2026 Tim Hanel
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

#include <core/logging.hpp>
#include <core/types.hpp>
#include <parse/CmdOptionsMap.hpp>
#include <parse/parser.hpp>
#include <random/random.hpp>
#include <utils/interpolation.hpp> /* interpolateWavelength*/

#include <algorithm>
#include <any>
#include <cassert>
#include <cmath> /* M_PI */
#include <cstdlib> /* exit() */
#include <filesystem>
#include <functional> /* bind, placeholders */
#include <limits>
#include <optional>
#include <ranges>
#include <sstream> /* stringstream */
#include <string> /* string */
#include <vector> /* vector */

namespace hase::parse
{
    namespace fs = std::filesystem;

    inline std::pair<std::string_view, std::optional<std::string_view>> splitEq(std::string_view arg)
    {
        auto const pos = arg.find('=');

        if(pos == std::string_view::npos)
        {
            return {arg, std::nullopt};
        }

        return {arg.substr(0, pos), arg.substr(pos + 1)};
    }

    struct WavelengthData
    {
        std::vector<double> sigmaAInterpolated;
        std::vector<double> sigmaEInterpolated;
        double maxSigmaA = 0;
        double maxSigmaE = 0;
    };

    /**
     *
     */
    hase::core::HostMesh createHostMeshFromFile(fs::path const rootPath)
    {
        // Parse experimentdata from files
        std::vector<unsigned> triangleNormalPoint = fileToVector<unsigned>(rootPath / "triangleNormalPoint.txt");
        std::vector<double> betaVolume = fileToVector<double>(rootPath / "betaVolume.txt");
        std::vector<int> forbiddenEdge = fileToVector<int>(rootPath / "forbiddenEdge.txt");
        std::vector<int> triangleNeighbors = fileToVector<int>(rootPath / "triangleNeighbors.txt");
        std::vector<double> triangleNormalsX = fileToVector<double>(rootPath / "triangleNormalsX.txt");
        std::vector<double> triangleNormalsY = fileToVector<double>(rootPath / "triangleNormalsY.txt");
        std::vector<double> triangleCenterX = fileToVector<double>(rootPath / "triangleCenterX.txt");
        std::vector<double> triangleCenterY = fileToVector<double>(rootPath / "triangleCenterY.txt");
        std::vector<double> points = fileToVector<double>(rootPath / "points.txt");
        std::vector<unsigned> trianglePointIndices = fileToVector<unsigned>(rootPath / "trianglePointIndices.txt");
        std::vector<float> triangleSurfaces = fileToVector<float>(rootPath / "triangleSurfaces.txt");
        unsigned numberOfPoints = fileToValue<unsigned>(rootPath / "numberOfPoints.txt");
        unsigned numberOfTriangles = fileToValue<unsigned>(rootPath / "numberOfTriangles.txt");
        unsigned numberOfLevels = fileToValue<unsigned>(rootPath / "numberOfLevels.txt");
        float thickness = fileToValue<float>(rootPath / "thickness.txt");
        float nTot = fileToValue<float>(rootPath / "nTot.txt");
        float crystalTFluo = fileToValue<float>(rootPath / "crystalTFluo.txt");
        unsigned claddingNumber = fileToValue<unsigned>(rootPath / "claddingNumber.txt");
        double claddingAbsorption = fileToValue<double>(rootPath / "claddingAbsorption.txt");
        std::vector<double> betaCells = fileToVector<double>(rootPath / "betaCells.txt");
        std::vector<unsigned> claddingCellTypes = fileToVector<unsigned>(rootPath / "claddingCellTypes.txt");
        std::vector<float> refractiveIndices = fileToVector<float>(rootPath / "refractiveIndices.txt");
        std::vector<float> reflectivities = fileToVector<float>(rootPath / "reflectivities.txt");
        return hase::core::HostMesh(
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
            claddingAbsorption);
    }

    bool validateHostMesh(hase::core::HostMesh const& mesh)
    {
        // assert input sizes
        assert(mesh.numberOfPoints == (mesh.points.size() / 2));
        assert(mesh.numberOfTriangles == mesh.trianglePointIndices.size() / 3);
        assert(mesh.triangleNormalPoint.size() == mesh.numberOfTriangles * 3);
        assert(mesh.triangleCenterY.size() == mesh.numberOfTriangles);
        assert(mesh.triangleCenterX.size() == mesh.numberOfTriangles);
        assert(mesh.triangleSurfaces.size() == mesh.numberOfTriangles);
        assert(mesh.betaVolume.size() == mesh.numberOfTriangles * (mesh.numberOfLevels - 1));
        assert(mesh.triangleNormalsX.size() == mesh.numberOfTriangles * 3);
        assert(mesh.triangleNormalsY.size() == mesh.numberOfTriangles * 3);
        assert(mesh.trianglePointIndices.size() == mesh.numberOfTriangles * 3);
        assert(mesh.forbiddenEdge.size() == mesh.numberOfTriangles * 3);
        assert(mesh.triangleNeighbors.size() == mesh.numberOfTriangles * 3);
        assert(mesh.betaCells.size() == mesh.numberOfPoints * mesh.numberOfLevels);
        assert(mesh.claddingCellTypes.size() == mesh.numberOfTriangles);
        assert(mesh.refractiveIndices.size() == 4);
        assert(mesh.reflectivities.size() == (mesh.refractiveIndices.size() / 2) * mesh.numberOfTriangles);
        assert(mesh.claddingCellTypes.size() == mesh.numberOfTriangles);

        // assert input data validity
        hase::core::assertRange(mesh.triangleNormalPoint, 0u, unsigned(mesh.numberOfPoints - 1), true);
        hase::core::assertMin(mesh.betaVolume, 0, false);
        hase::core::assertRange(mesh.forbiddenEdge, -1, 2, true);
        hase::core::assertRange(mesh.triangleNeighbors, -1, int(mesh.numberOfTriangles - 1), true);
        hase::core::assertRange(mesh.triangleNormalsX, -1, 1, false);
        hase::core::assertRange(mesh.triangleNormalsY, -1, 1, false);

        hase::core::assertRange(
            mesh.triangleCenterX,
            *std::min_element(mesh.points.begin(), mesh.points.end()),
            *std::max_element(mesh.points.begin(), mesh.points.end()),
            false);

        hase::core::assertRange(
            mesh.triangleCenterY,
            *std::min_element(mesh.points.begin(), mesh.points.end()),
            *std::max_element(mesh.points.begin(), mesh.points.end()),
            false);

        hase::core::assertRange(mesh.trianglePointIndices, 0u, unsigned(mesh.numberOfPoints - 1), true);
        hase::core::assertMin(mesh.triangleSurfaces, 0, false);
        hase::core::assertMin(mesh.betaCells, 0, false);
        hase::core::assertRange(mesh.refractiveIndices, 0, 5, false);
        hase::core::assertRange(mesh.reflectivities, 0, 1, false);

        return true;
    }

    void validateSampleRange(hase::core::ComputeParameters& compute, unsigned const numberOfSamples)
    {
        bool const minSet = compute.minSampleRange != std::numeric_limits<unsigned>::max();
        bool const maxSet = compute.maxSampleRange != std::numeric_limits<unsigned>::max();

        if(!minSet && !maxSet)
        {
            hase::core::dout(V_INFO) << hase::core::CompSwitch::min_sample_i << "/"
                                     << hase::core::CompSwitch::max_sample_i
                                     << " not set! Assuming a sampling point range of 0 to " << (numberOfSamples - 1)
                                     << std::endl;
            compute.minSampleRange = 0;
            compute.maxSampleRange = numberOfSamples - 1;
            return;
        }

        if(minSet != maxSet)
        {
            throw std::runtime_error(
                "Options " + std::string(hase::core::CompSwitch::min_sample_i) + "/"
                + std::string(hase::core::CompSwitch::max_sample_i) + " must be used together!");
        }

        if(compute.maxSampleRange < compute.minSampleRange)
        {
            throw std::runtime_error(
                std::string(hase::core::CompSwitch::max_sample_i) + " < "
                + std::string(hase::core::CompSwitch::min_sample_i) + "!");
        }

        if(compute.minSampleRange >= numberOfSamples)
        {
            throw std::runtime_error(
                std::string(hase::core::CompSwitch::min_sample_i) + " is out of range! (There are only "
                + std::to_string(numberOfSamples) + " sampling points)");
        }

        if(compute.maxSampleRange >= numberOfSamples)
        {
            throw std::runtime_error(
                std::string(hase::core::CompSwitch::max_sample_i) + " is out of range! (There are only "
                + std::to_string(numberOfSamples) + " sampling points)");
        }

        unsigned const samplesForNode = compute.maxSampleRange - compute.minSampleRange + 1;
        if(compute.numDevices > samplesForNode)
        {
            hase::core::dout(V_WARNING) << "More GPUs requested than there are sample points. "
                                        << "Number of used GPUs reduced to " << samplesForNode << std::endl;
            compute.numDevices = samplesForNode;
        }
        hase::core::dout(V_WARNING)
            << "Manual sampling range set via " << hase::core::CompSwitch::min_sample_i << "/"
            << hase::core::CompSwitch::max_sample_i
            << ". This may allocate more memory than actually required. "
               "Normally, the sampling range is determined automatically from the number of mesh "
               "points and the number of layers, so setting it manually is not required for most "
               "use-cases."
            << std::endl;
    }

    inline void ensureOrThrow(bool cond, std::string const& msg)
    {
        if(!cond)
            throw std::runtime_error(msg);
    }

    void validateInput(hase::core::ExperimentParameters& experiment, hase::core::ComputeParameters& compute)
    {
        // if(compute.backend != hase::core::Backend::CPU && compute.backend != hase::core::Backend::GPU)
        // {
        //     throw std::runtime_error(
        //         hase::core::CompSwitch::backend + " must be either \"" +
        //         std::string(hase::core::Backend::CPU) + "\" or \"" +
        //         std::string(hase::core::Backend::GPU) + "\"");
        // }

        if(compute.backend == hase::core::Backend::CPU && compute.parallelMode != hase::core::ParallelMode::SINGLE)
        {
            hase::core::dout(V_WARNING) << hase::core::CompSwitch::backend << " \"" << hase::core::Backend::CPU
                                        << "\" does only support " << hase::core::CompSwitch::parallel_mode
                                        << " \"single\"! (will be ignored)" << std::endl;
        }

        if(compute.parallelMode != hase::core::ParallelMode::SINGLE
#if !defined(DISABLE_MPI)
           && compute.parallelMode != hase::core::ParallelMode::MPI
#endif
        )
        {
            std::stringstream ss;
            ss << hase::core::CompSwitch::parallel_mode << " must be one of \"" << hase::core::ParallelMode::SINGLE
               << "\"";
#if !defined(DISABLE_MPI)
            ss << ", \"" << hase::core::ParallelMode::MPI << "\"";
#endif
            throw std::runtime_error(ss.str());
        }

        ensureOrThrow(
            experiment.minRaysPerSample > 0,
            "Please specify a positive value for --" + hase::core::ExpSwitch::min_rays);

        if(experiment.maxRaysPerSample < experiment.minRaysPerSample)
        {
            hase::core::dout(V_WARNING) << hase::core::ExpSwitch::max_rays << " < " << hase::core::ExpSwitch::min_rays
                                        << ". Auto-increasing " << hase::core::ExpSwitch::max_rays
                                        << " (will be non-adaptive!)" << std::endl;
            experiment.maxRaysPerSample = experiment.minRaysPerSample;
        }


        if(hase::core::verbosity >= 64)
        {
            hase::core::verbosity = 63;
            hase::core::dout(V_WARNING) << "Verbosity level should be between 0 (quiet) and 63 (all). "
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
            ensureOrThrow(
                fs::exists(compute.inputPath) && fs::is_directory(compute.inputPath),
                "The specified input path does not exist, is no directory, or has insufficient permissions.");

            ensureOrThrow(!fs::is_empty(compute.inputPath), "The specified input folder is empty.");
        }

        if(!compute.outputPath.empty())
        {
            ensureOrThrow(
                fs::exists(compute.outputPath) && fs::is_directory(compute.outputPath),
                "The specified output path does not exist (or permission denied).");
        }
    }

    WavelengthData calculateSigmas(
        std::vector<double>& sigmaA,
        std::vector<double>& sigmaE,
        std::vector<double>& lambdaA,
        std::vector<double>& lambdaE,
        unsigned lambdaResolution,
        bool monochromatic = false)
    {
        WavelengthData waveD;

        if(monochromatic)
        {
            if(sigmaA.empty() || sigmaE.empty())
            {
                throw std::runtime_error("Cannot use monochromatic propagation without sigmaA/sigmaE values.");
            }

            waveD.sigmaAInterpolated = {sigmaA.front()};
            waveD.sigmaEInterpolated = {sigmaE.front()};
            waveD.maxSigmaA = sigmaA.front();
            waveD.maxSigmaE = sigmaE.front();
            return waveD;
        }

        assert(sigmaA.size() == lambdaA.size());
        assert(sigmaE.size() == lambdaE.size());

        bool const singleValueInput
            = sigmaA.size() == 1 && sigmaE.size() == 1 && lambdaA.size() == 1 && lambdaE.size() == 1;
        if(singleValueInput)
        {
            waveD.sigmaAInterpolated = {sigmaA.front()};
            waveD.sigmaEInterpolated = {sigmaE.front()};
            waveD.maxSigmaA = sigmaA.front();
            waveD.maxSigmaE = sigmaE.front();
            return waveD;
        }

        lambdaResolution = std::max(lambdaResolution, (unsigned) lambdaA.size());
        lambdaResolution = std::max(lambdaResolution, (unsigned) lambdaE.size());

        // Interpolate sigmaA / sigmaE function
        waveD.sigmaAInterpolated = hase::utils::interpolateLinear(sigmaA, lambdaA, lambdaResolution);
        waveD.sigmaEInterpolated = hase::utils::interpolateLinear(sigmaE, lambdaE, lambdaResolution);
        assert(waveD.sigmaAInterpolated.size() == waveD.sigmaEInterpolated.size());

        // Calc max sigmaA / sigmaE
        for(unsigned i = 0; i < sigmaE.size(); ++i)
        {
            if(sigmaE.at(i) > waveD.maxSigmaE)
            {
                waveD.maxSigmaE = sigmaE.at(i);
                waveD.maxSigmaA = sigmaA.at(i);
            }
        }
        return waveD;
    }

    WavelengthData calculateSigmas(fs::path inputPath, unsigned lambdaResolution, bool monochromatic)
    {
        std::vector<double> sigmaA = fileToVector<double>(inputPath / "sigmaA.txt");
        std::vector<double> sigmaE = fileToVector<double>(inputPath / "sigmaE.txt");

        if(monochromatic)
        {
            std::vector<double> lambdaA;
            std::vector<double> lambdaE;
            return calculateSigmas(sigmaA, sigmaE, lambdaA, lambdaE, lambdaResolution, monochromatic);
        }

        // Parse wavelengths from files
        std::vector<double> lambdaA = fileToVector<double>(inputPath / "lambdaA.txt");
        std::vector<double> lambdaE = fileToVector<double>(inputPath / "lambdaE.txt");
        return calculateSigmas(sigmaA, sigmaE, lambdaA, lambdaE, lambdaResolution);
    }

    void checkPositive(int i, std::string const& name, int threshhold = 0)
    {
        if(i < threshhold)
        {
            hase::core::verbosity |= V_ERROR;
            hase::core::dout(V_ERROR) << name << " must have a positive argument!" << std::endl;
            throw std::invalid_argument(std::to_string(i));
        }
    }

    CmdOptionsMap::Option const& selectedNumDevicesOption(CmdOptionsMap const& vm)
    {
        auto const numDevicesSet = vm.isExplicitlySet(hase::core::CompSwitch::numDevices);

        auto const explicitCount = static_cast<int>(numDevicesSet);
        if(explicitCount > 1)
        {
            throw std::runtime_error("Use only --" + hase::core::CompSwitch::numDevices);
        }
        return vm.options.at(hase::core::CompSwitch::numDevices);
    }

    void initializeResult(hase::core::Result& result, unsigned const numberOfSamples)
    {
        result = hase::core::Result(
            std::vector<float>(numberOfSamples, 0.0f),
            std::vector<double>(numberOfSamples, std::numeric_limits<double>::max()),
            std::vector<unsigned>(numberOfSamples, 0u),
            std::vector<double>(numberOfSamples, 0.0),
            std::vector<unsigned>(numberOfSamples, 0u));
    }

    int pythonParse(
        hase::core::ExperimentParameters& experiment,
        hase::core::ComputeParameters& compute,
        hase::core::HostMesh& host_mesh,
        hase::core::Result& result)
    {
        host_mesh.calcTotalReflectionAngles();
        validateHostMesh(host_mesh);
        validateInput(experiment, compute);
        unsigned numberOfSamples = host_mesh.numberOfPoints * host_mesh.numberOfLevels;

        validateSampleRange(compute, numberOfSamples);
        auto wave = calculateSigmas(
            experiment.sigmaA,
            experiment.sigmaE,
            experiment.lambdaA,
            experiment.lambdaE,
            experiment.spectral,
            experiment.monochromatic);
        experiment.sigmaA = wave.sigmaAInterpolated;
        experiment.sigmaE = wave.sigmaEInterpolated;
        experiment.maxSigmaA = wave.maxSigmaA;
        experiment.maxSigmaE = wave.maxSigmaE;

        initializeResult(result, numberOfSamples);

        return 0;
    }

    inline CmdOptionsMap constructOptionsMapFromFile(std::filesystem::path const& cfgPath)
    {
        CmdOptionsMap cfg{};

        std::ifstream file(cfgPath);
        if(!file)
        {
            throw std::runtime_error("Could not open config file: " + cfgPath.string());
        }

        std::string line;
        std::size_t lineNo = 0;

        while(std::getline(file, line))
        {
            ++lineNo;

            auto trim = [](std::string_view s) -> std::string_view
            {
                auto isWs = [](unsigned char c) { return std::isspace(c); };

                while(!s.empty() && isWs(static_cast<unsigned char>(s.front())))
                    s.remove_prefix(1);

                while(!s.empty() && isWs(static_cast<unsigned char>(s.back())))
                    s.remove_suffix(1);

                return s;
            };

            std::string_view view = trim(line);

            if(view.empty() || view.starts_with("#") || view.starts_with(";"))
                continue;

            auto pos = view.find('=');
            if(pos == std::string_view::npos)
            {
                throw std::runtime_error(
                    "Invalid config line " + std::to_string(lineNo) + " in '" + cfgPath.string()
                    + "': expected key=value");
            }

            auto key = trim(view.substr(0, pos));
            auto value = trim(view.substr(pos + 1));

            if(key.empty())
            {
                throw std::runtime_error(
                    "Invalid config line " + std::to_string(lineNo) + " in '" + cfgPath.string() + "': empty key");
            }

            cfg.set(std::string(key), std::string(value));
        }

        cfg.checkRequired();
        return cfg;
    }

    int constructSimulationContextFromOptionsMap(
        CmdOptionsMap& vm,
        hase::core::ExperimentParameters& experiment,
        hase::core::ComputeParameters& compute,
        std::optional<hase::core::HostMesh>& hostMesh,
        hase::core::Result& result)
    {
        // Set/Test device to run experiment with
        // TODO: this call takes a LOT of time (2-5s). Can this be avoided?
        //      maybe move this to a place where GPUs are actually needed (for_loops_clad doesn't even need GPUs!)
        auto const& nr_gpuOpt = selectedNumDevicesOption(vm);
        checkPositive(nr_gpuOpt.as<int>(), nr_gpuOpt.name);
        hostMesh.emplace(createHostMeshFromFile(vm[hase::core::ExpSwitch::input_path].as<fs::path>()));
        hostMesh->calcTotalReflectionAngles();
        validateHostMesh(*hostMesh);
        bool const monochromatic = vm[hase::core::ExpSwitch::monochromatic].as<bool>();
        WavelengthData waveD = calculateSigmas(
            vm[hase::core::ExpSwitch::input_path].as<fs::path>(),
            vm.isExplicitlySet(hase::core::ExpSwitch::spectral) ? vm[hase::core::ExpSwitch::spectral].as<int>() : 0,
            monochromatic);

        auto& min_rayOpt = vm[hase::core::ExpSwitch::min_rays];
        auto& max_rayOpt = vm[hase::core::ExpSwitch::max_rays];
        checkPositive(min_rayOpt.as<int>(), min_rayOpt.name);
        checkPositive(max_rayOpt.as<int>(), max_rayOpt.name);
        experiment = hase::core::ExperimentParameters(
            static_cast<unsigned>(min_rayOpt.as<int>()),
            static_cast<unsigned>(max_rayOpt.as<int>()),
            waveD.sigmaAInterpolated,
            waveD.sigmaEInterpolated,
            waveD.maxSigmaA,
            waveD.maxSigmaE,
            vm[hase::core::ExpSwitch::mse].as<double>(),
            vm[hase::core::ExpSwitch::reflection].as<bool>(),
            monochromatic);

        auto& repetionOpt = vm[hase::core::CompSwitch::repetitions];
        auto& adaptiveStepOpt = vm[hase::core::CompSwitch::adaptive_steps];
        checkPositive(repetionOpt.as<int>(), repetionOpt.name);
        checkPositive(adaptiveStepOpt.as<int>(), adaptiveStepOpt.name);
        compute = hase::core::ComputeParameters(
            static_cast<unsigned>(repetionOpt.as<int>()),
            static_cast<unsigned>(adaptiveStepOpt.as<int>()),
            0,
            vm[hase::core::CompSwitch::backend].as<std::string>(),
            vm[hase::core::CompSwitch::parallel_mode].as<std::string>(),
            vm[hase::core::CompSwitch::write_vtk].as<bool>(),
            vm[hase::core::ExpSwitch::input_path].as<fs::path>(),
            vm[hase::core::ExpSwitch::output_path].as<fs::path>(),
            std::vector<unsigned>{},
            vm[hase::core::CompSwitch::min_sample_i].as<unsigned>(),
            vm[hase::core::CompSwitch::max_sample_i].as<unsigned>(),
            nr_gpuOpt.as<int>());
        if(vm.isExplicitlySet(hase::core::CompSwitch::rng_seed))
        {
            compute.rngSeed = vm[hase::core::CompSwitch::rng_seed].as<unsigned>();
        }

        validateInput(experiment, compute);
        validateSampleRange(compute, hostMesh->numberOfLevels * hostMesh->numberOfPoints);
        initializeResult(result, hostMesh->numberOfLevels * hostMesh->numberOfPoints);
        return 0;
    }

    int parseFile(
        std::filesystem::path const& file,
        hase::core::ExperimentParameters& experiment,
        hase::core::ComputeParameters& compute,
        std::optional<hase::core::HostMesh>& hostMesh,
        hase::core::Result& result)
    {
        auto vm = constructOptionsMapFromFile(file);
        constructSimulationContextFromOptionsMap(vm, experiment, compute, hostMesh, result);
        return 0;
    }

    int parse(
        int const argc,
        char** argv,
        hase::core::ExperimentParameters& experiment,
        hase::core::ComputeParameters& compute,
        std::optional<hase::core::HostMesh>& hostMesh,
        hase::core::Result& result)
    {
        for(int i = 1; i < argc; ++i)
        {
            std::string_view arg = argv[i];
            auto const [name, inlineValue] = splitEq(arg);

            if(name == "-cfg" || name == "--cfg")
            {
                std::filesystem::path cfgPath;

                if(inlineValue.has_value())
                {
                    cfgPath = std::filesystem::path{std::string(*inlineValue)};
                }
                else
                {
                    if(i + 1 >= argc)
                    {
                        throw std::runtime_error("Missing value for option '" + std::string(name) + "'");
                    }

                    cfgPath = std::filesystem::path{argv[i + 1]};
                }

                hase::core::dout(V_INFO) << "Parsing simulation options using .cfg file '" << cfgPath.string()
                                         << "' -- other command line options are ignored!" << std::endl;

                parseFile(cfgPath, experiment, compute, hostMesh, result);
                return 0;
            }
        }

        auto vm = parseCommandLine(argc, argv);
        printCommandLine(vm);
        constructSimulationContextFromOptionsMap(vm, experiment, compute, hostMesh, result);

        return 0;
    }

    inline CmdOptionsMap parseCommandLine(int argc, char** argv)
    {
        CmdOptionsMap cfg{};

        for(int i = 1; i < argc; ++i)
        {
            std::string_view arg = argv[i];

            // force --name=value form
            auto splitArg = splitEq(arg);
            auto name = splitArg.first;
            auto inlineValue = splitArg.second;

            auto getValue = [&](std::string_view optName) -> std::string_view
            {
                if(inlineValue.has_value())
                {
                    return *inlineValue;
                }
                return requireValue(i, argc, argv, optName);
            };
            auto removePreDashes = [&](std::string_view optName) -> std::string_view
            {
                if(optName.empty())
                {
                    throw std::runtime_error("Invalid option name: empty string");
                }

                if(optName[0] != '-')
                {
                    throw std::runtime_error(
                        "Invalid option name '" + std::string(optName) + "': expected leading '-' or '--'");
                }

                auto const pos = optName.find_first_not_of('-');
                if(pos == std::string_view::npos)
                {
                    throw std::runtime_error(
                        "Invalid option name '" + std::string(optName) + "': missing option name after dashes");
                }

                return optName.substr(pos);
            };
            auto cleanName = removePreDashes(name);
            if(cleanName == "help" || cleanName == "h")
            {
                cfg.printDescription();
                std::terminate();
            }
            cfg.set(std::string(cleanName), std::string(getValue(name)));
        }
        cfg.checkRequired();
        return cfg;
    }

    void printCommandLine(CmdOptionsMap const& vm)
    {
        for(auto const& [name, opt] : vm.options)
        {
            std::string const& shownValue = opt.hasValue ? opt.value : opt.defaultValue;
            hase::core::dout(V_INFO) << name << ": " << shownValue << std::endl;
        }
    }

} // namespace hase::parse
