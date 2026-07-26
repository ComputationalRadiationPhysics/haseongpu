#include <catch2/catch_test_macros.hpp>
#include <openpmd/SimulationSnapshotWriter.hpp>

#include <thread>
#include <vector>

TEST_CASE("simulation snapshot writer runs synchronously when requested", "[openpmd][mpi]")
{
    auto const callerThread = std::this_thread::get_id();
    std::thread::id writerThread;
    std::vector<unsigned> writtenSteps;
    hase::openpmd::AsyncSimulationSnapshotWriter writer{
        true,
        [&](hase::core::SimulationSnapshot const& snapshot)
        {
            writerThread = std::this_thread::get_id();
            writtenSteps.push_back(snapshot.step);
        },
        false};

    hase::core::SimulationSnapshot snapshot;
    snapshot.step = 7u;
    writer.enqueue(snapshot);

    REQUIRE(writtenSteps == std::vector<unsigned>{7u});
    CHECK(writerThread == callerThread);
    CHECK_NOTHROW(writer.finish());
}

TEST_CASE("simulation snapshot writer retains asynchronous mode", "[openpmd]")
{
    auto const callerThread = std::this_thread::get_id();
    std::thread::id writerThread;
    std::vector<unsigned> writtenSteps;
    hase::openpmd::AsyncSimulationSnapshotWriter writer{
        true,
        [&](hase::core::SimulationSnapshot const& snapshot)
        {
            writerThread = std::this_thread::get_id();
            writtenSteps.push_back(snapshot.step);
        }};

    hase::core::SimulationSnapshot first;
    first.step = 3u;
    hase::core::SimulationSnapshot second;
    second.step = 4u;
    writer.enqueue(first);
    writer.enqueue(second);
    writer.finish();

    REQUIRE(writtenSteps == std::vector<unsigned>{3u, 4u});
    CHECK(writerThread != callerThread);
}
