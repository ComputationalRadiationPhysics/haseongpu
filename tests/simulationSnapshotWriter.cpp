#include <catch2/catch_test_macros.hpp>
#include <core/simulationRunControl.hpp>
#include <openpmd/SimulationSnapshotWriter.hpp>

#include <chrono>
#include <future>
#include <thread>
#include <type_traits>
#include <vector>

namespace
{
    template<typename T>
    concept ContainsMeshMember = requires(T value) { value.mesh; };
    template<typename T>
    concept ContainsPointBetaMember = requires(T value) { value.betaCells; };
} // namespace

TEST_CASE("simulation snapshots contain dynamic state instead of a host mesh", "[openpmd]")
{
    STATIC_REQUIRE_FALSE(ContainsMeshMember<hase::core::SimulationSnapshot>);
    STATIC_REQUIRE_FALSE(ContainsPointBetaMember<hase::core::SimulationSnapshot>);
    STATIC_REQUIRE(std::is_same_v<decltype(hase::core::SimulationSnapshot::betaVolume), std::vector<double>>);
    STATIC_REQUIRE(std::is_same_v<decltype(hase::core::SimulationSnapshot::aseResult), hase::core::Result>);
}

TEST_CASE("simulation run fields expose only cell-centered state", "[simulation]")
{
    using hase::core::SimulationControlField;
    using hase::core::SimulationOutputField;

    CHECK(
        SimulationOutputField::all()
        == std::vector<std::string>{
            "beta_volume",
            "phi_ase",
            "standard_error",
            "relative_standard_error",
            "total_rays",
            "dndt_ase",
            "dndt_pump"});
    CHECK(SimulationControlField::all() == std::vector<std::string>{"beta_volume"});
}

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

TEST_CASE("simulation snapshot writer bounds pending output", "[openpmd]")
{
    using namespace std::chrono_literals;

    std::promise<void> firstWriteStarted;
    auto firstWriteStartedFuture = firstWriteStarted.get_future();
    std::promise<void> releaseFirstWrite;
    auto releaseFirstWriteFuture = releaseFirstWrite.get_future().share();
    std::vector<unsigned> writtenSteps;
    hase::openpmd::AsyncSimulationSnapshotWriter writer{
        true,
        [&](hase::core::SimulationSnapshot const& snapshot)
        {
            if(snapshot.step == 1u)
            {
                firstWriteStarted.set_value();
                releaseFirstWriteFuture.wait();
            }
            writtenSteps.push_back(snapshot.step);
        }};

    hase::core::SimulationSnapshot first;
    first.step = 1u;
    hase::core::SimulationSnapshot second;
    second.step = 2u;
    hase::core::SimulationSnapshot third;
    third.step = 3u;

    writer.enqueue(first);
    auto const writerStarted = firstWriteStartedFuture.wait_for(2s) == std::future_status::ready;
    CHECK(writerStarted);
    if(!writerStarted)
    {
        releaseFirstWrite.set_value();
        writer.finish();
        return;
    }

    writer.enqueue(second);
    auto thirdEnqueue = std::async(std::launch::async, [&] { writer.enqueue(third); });
    auto const blockedWithOnePending = thirdEnqueue.wait_for(100ms) == std::future_status::timeout;
    releaseFirstWrite.set_value();
    auto const completedAfterDrain = thirdEnqueue.wait_for(2s) == std::future_status::ready;

    CHECK(blockedWithOnePending);
    CHECK(completedAfterDrain);
    CHECK_NOTHROW(thirdEnqueue.get());
    CHECK_NOTHROW(writer.finish());
    CHECK(writtenSteps == std::vector<unsigned>{1u, 2u, 3u});
}
