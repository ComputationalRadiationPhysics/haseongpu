#pragma once

#include <core/mesh.hpp>
#include <core/types.hpp>

#include <cstdint>
#include <filesystem>
#include <functional>
#include <string>

#if defined(MPI_FOUND) && !defined(DISABLE_MPI)
#    include <mpi.h>
#endif

namespace openPMD
{
    class Iteration;
    class Series;
} // namespace openPMD

namespace hase::openpmd
{

    class Parser
    {
    public:
        Parser(std::filesystem::path inputPath, std::filesystem::path outputPath);

#if defined(MPI_FOUND) && !defined(DISABLE_MPI)
        Parser(std::filesystem::path inputPath, std::filesystem::path outputPath, MPI_Comm comm);
#endif

        [[nodiscard]] core::SimulationContext read();
        void writeResult(core::Result const& result, core::HostMesh const& mesh);
        void processAll(std::function<void(core::SimulationContext&)> process);

    private:
        [[nodiscard]] bool isHeadRank() const;
        [[nodiscard]] bool hasStaticMeshUpdate(openPMD::Iteration const& iteration) const;
        [[nodiscard]] bool containsStaticMeshUpdate(openPMD::Iteration const& iteration) const;
        void validateDynamicOnlyIteration(
            openPMD::Iteration const& iteration,
            core::SimulationContext const& simulation) const;
        [[nodiscard]] core::SimulationContext readIteration(openPMD::Series& series, openPMD::Iteration& iteration);
        void updateDynamicIteration(
            openPMD::Series& series,
            openPMD::Iteration& iteration,
            core::SimulationContext& simulation);
        void writeResultIteration(
            openPMD::Series& series,
            std::uint64_t iterationIndex,
            core::Result const& result,
            core::HostMesh const& mesh);

        std::filesystem::path m_inputPath;
        std::filesystem::path m_outputPath;
        std::string m_meshGroup = "core";

#if defined(MPI_FOUND) && !defined(DISABLE_MPI)
        MPI_Comm m_comm = MPI_COMM_WORLD;
#endif
    };

} // namespace hase::openpmd
