#include <core/simulation.hpp>
#include <openpmd/OpenPmdParser.hpp>

#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>

namespace
{
    struct Paths
    {
        std::filesystem::path input;
        std::filesystem::path output;
        bool runSimulation = false;
    };

    std::optional<std::string_view> valueFor(std::string_view arg, std::string_view name)
    {
        std::string const prefix = "--" + std::string(name) + "=";
        if(arg.starts_with(prefix))
        {
            return arg.substr(prefix.size());
        }
        return std::nullopt;
    }

    Paths parsePaths(int argc, char** argv)
    {
        Paths paths;
        for(int i = 1; i < argc; ++i)
        {
            std::string_view const arg = argv[i];
            if(auto value = valueFor(arg, "input-path"))
            {
                paths.input = std::string(*value);
                continue;
            }
            if(auto value = valueFor(arg, "output-path"))
            {
                paths.output = std::string(*value);
                continue;
            }
            if(arg == "--run-simulation")
            {
                paths.runSimulation = true;
                continue;
            }
            throw std::runtime_error(
                "Unsupported argument '" + std::string(arg)
                + "'. calcPhiASE only accepts --input-path, --output-path, and --run-simulation.");
        }

        if(paths.input.empty())
        {
            throw std::runtime_error("Missing required --input-path=<openPMD-series>.");
        }
        if(paths.output.empty())
        {
            throw std::runtime_error("Missing required --output-path=<openPMD-series>.");
        }
        return paths;
    }
} // namespace

int main(int argc, char** argv)
{
    try
    {
        auto paths = parsePaths(argc, argv);
#if defined(MPI_FOUND) && !defined(DISABLE_MPI)
        MPI_Init(&argc, &argv);
        hase::openpmd::Parser openPmdParser{paths.input, paths.output, MPI_COMM_WORLD};
#else
        hase::openpmd::Parser openPmdParser{paths.input, paths.output};
#endif

        if(paths.runSimulation)
        {
            openPmdParser.runTimeSteppedSimulation();
        }
        else
        {
            openPmdParser.processAll(
                [](hase::core::SimulationContext& simulation)
                {
                    int const result = hase::core::startSimulation<false>(
                        simulation.experiment,
                        simulation.compute,
                        simulation.result,
                        simulation.mesh);
                    if(result != 0)
                    {
                        throw std::runtime_error("simulation failed with return code " + std::to_string(result));
                    }
                });
        }

#if defined(MPI_FOUND) && !defined(DISABLE_MPI)
        MPI_Finalize();
#endif
        return 0;
    }
    catch(std::exception const& e)
    {
        std::cerr << "calcPhiASE failed: " << e.what() << '\n';
#if defined(MPI_FOUND) && !defined(DISABLE_MPI)
        int mpiInitialized = 0;
        MPI_Initialized(&mpiInitialized);
        int mpiFinalized = 0;
        MPI_Finalized(&mpiFinalized);
        if(mpiInitialized && !mpiFinalized)
        {
            MPI_Finalize();
        }
#endif
        return 1;
    }
}
