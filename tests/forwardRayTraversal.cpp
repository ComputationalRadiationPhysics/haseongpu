#include <alpaka/math.hpp>

#include <alpakaUtils/DevBundle.hpp>
#include <alpakaUtils/memory.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <core/mesh.hpp>
#include <kernels/forward/rayTransition.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <ranges>
#include <type_traits>
#include <vector>

namespace hase::tests
{
    using TestBackends = std::decay_t<
        decltype(alpaka::onHost::allBackends(alpaka::onHost::enabledApis, alpaka::exec::enabledExecutors))>;

    hase::core::HostMesh makeTraversalMesh(
        std::vector<hase::core::Point> const& points,
        std::vector<std::array<unsigned, hase::core::tet4VertexCount>> const& cells)
    {
        hase::core::HostMesh mesh;
        mesh.numberOfCells = static_cast<unsigned>(cells.size());
        mesh.numberOfMeshPoints = static_cast<unsigned>(points.size());
        mesh.numberOfCellVertices = hase::core::tet4VertexCount;
        mesh.numberOfFacesPerCell = hase::core::tet4FaceCount;
        mesh.points.resize(points.size() * 3u);
        for(unsigned point = 0u; point < points.size(); ++point)
        {
            mesh.points[point] = points[point].x;
            mesh.points[point + points.size()] = points[point].y;
            mesh.points[point + 2u * points.size()] = points[point].z;
        }

        using Face = std::array<unsigned, hase::core::tet4FaceWidth>;
        std::vector<Face> faces;
        faces.reserve(cells.size() * hase::core::tet4FaceCount);
        for(auto const& cell : cells)
        {
            mesh.cellPointIndices.insert(mesh.cellPointIndices.end(), cell.cbegin(), cell.cend());
            for(unsigned localFace = 0u; localFace < hase::core::tet4FaceCount; ++localFace)
            {
                Face face{};
                unsigned faceVertex = 0u;
                for(unsigned localVertex = 0u; localVertex < hase::core::tet4VertexCount; ++localVertex)
                {
                    if(localVertex != localFace)
                    {
                        face[faceVertex++] = cell[localVertex];
                    }
                }
                std::ranges::sort(face);
                faces.push_back(face);
                mesh.cellFaces.insert(mesh.cellFaces.end(), face.cbegin(), face.cend());
            }
        }

        mesh.cellNeighborCells.assign(faces.size(), -1);
        mesh.cellNeighborLocalFaces.assign(faces.size(), -1);
        for(unsigned face = 0u; face < faces.size(); ++face)
        {
            for(unsigned candidate = face + 1u; candidate < faces.size(); ++candidate)
            {
                if(faces[face] != faces[candidate])
                {
                    continue;
                }
                unsigned const cell = face / hase::core::tet4FaceCount;
                unsigned const localFace = face % hase::core::tet4FaceCount;
                unsigned const neighbor = candidate / hase::core::tet4FaceCount;
                unsigned const neighborLocalFace = candidate % hase::core::tet4FaceCount;
                mesh.cellNeighborCells[face] = static_cast<int>(neighbor);
                mesh.cellNeighborLocalFaces[face] = static_cast<int>(neighborLocalFace);
                mesh.cellNeighborCells[candidate] = static_cast<int>(cell);
                mesh.cellNeighborLocalFaces[candidate] = static_cast<int>(localFace);
            }
        }
        mesh.precomputeBarycentricFacePlanes();
        return mesh;
    }

    struct TraversalResult
    {
        hase::kernels::forward::Tet4FaceIntersection intersection;
        hase::core::Point intersectionPoint;
        hase::kernels::forward::Tet4FaceTransition transition;
        hase::kernels::forward::Tet4FaceIntersection nextIntersection;
    };

    static_assert(std::is_trivially_copyable_v<TraversalResult>);

    struct TraverseOneRay
    {
        template<typename TAcc, alpaka::concepts::IMdSpan TResult>
        ALPAKA_FN_ACC void operator()(
            TAcc const&,
            hase::core::DeviceMeshView const mesh,
            unsigned const startCell,
            hase::core::Point const origin,
            hase::core::Point const direction,
            TResult result) const
        {
            TraversalResult traversal;
            traversal.intersection
                = hase::kernels::forward::nextFaceIntersection(mesh, startCell, origin, direction, -1);
            traversal.intersectionPoint
                = hase::kernels::forward::advance(origin, direction, traversal.intersection.length);
            traversal.transition = hase::kernels::forward::transitionAcrossIntersection(
                mesh,
                startCell,
                traversal.intersection,
                traversal.intersectionPoint,
                direction);
            if(traversal.transition.status == hase::kernels::forward::Tet4TransitionStatus::enteredCell)
            {
                traversal.nextIntersection = hase::kernels::forward::nextFaceIntersection(
                    mesh,
                    traversal.transition.cell,
                    traversal.intersectionPoint,
                    direction,
                    traversal.transition.forbiddenFace);
            }
            result[0u] = traversal;
        }
    };

    TraversalResult traverseOneRay(
        auto& device,
        auto const executor,
        hase::core::HostMesh& hostMesh,
        unsigned const startCell,
        hase::core::Point const origin,
        hase::core::Point const direction)
    {
        auto queue = device.makeQueue(alpaka::queueKind::blocking);
        auto deviceMesh = hostMesh.toDevice(device);
        std::vector<TraversalResult> result(1u);
        auto deviceResult = hase::alpakaUtils::toDevice(queue, result);
        auto const frameSpec = hase::alpakaUtils::getFrameSpec<std::uint32_t>(device, executor, alpaka::Vec{1u});
        queue.enqueue(
            frameSpec,
            alpaka::KernelBundle{TraverseOneRay{}, deviceMesh.toView(), startCell, origin, direction, deviceResult});
        alpaka::onHost::memcpy(queue, result, deviceResult);
        alpaka::onHost::wait(queue);
        return result.front();
    }

} // namespace hase::tests

TEMPLATE_LIST_TEST_CASE(
    "forward ray crosses a vertex shared by 30 tetrahedra",
    "[forward][traversal][edge-case]",
    hase::tests::TestBackends)
{
    auto const backend = TestType::makeDict();
    auto deviceSelector = alpaka::onHost::makeDeviceSelector(backend[alpaka::object::deviceSpec]);
    if(!deviceSelector.isAvailable())
    {
        SUCCEED("No device available for " << backend[alpaka::object::deviceSpec].getName());
        return;
    }
    auto device = deviceSelector.makeDevice(0);

    constexpr unsigned tetrahedraPerHalf = 15u;
    constexpr unsigned sharedVertex = 0u;
    constexpr double pi = 3.14159265358979323846;
    std::vector<hase::core::Point> points{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}};
    for(unsigned point = 0u; point < tetrahedraPerHalf; ++point)
    {
        double const angle = 2.0 * pi * static_cast<double>(point) / tetrahedraPerHalf;
        points.push_back({std::cos(angle), std::sin(angle), 0.0});
    }

    std::vector<std::array<unsigned, hase::core::tet4VertexCount>> cells;
    cells.reserve(2u * tetrahedraPerHalf);
    for(unsigned cell = 0u; cell < tetrahedraPerHalf; ++cell)
    {
        unsigned const ring0 = 3u + cell;
        unsigned const ring1 = 3u + (cell + 1u) % tetrahedraPerHalf;
        cells.push_back({sharedVertex, 1u, ring0, ring1});
    }
    for(unsigned cell = 0u; cell < tetrahedraPerHalf; ++cell)
    {
        unsigned const ring0 = 3u + cell;
        unsigned const ring1 = 3u + (cell + 1u) % tetrahedraPerHalf;
        cells.push_back({sharedVertex, 2u, ring1, ring0});
    }
    auto mesh = hase::tests::makeTraversalMesh(points, cells);

    CHECK(mesh.numberOfCells == 30u);
    CHECK(
        std::ranges::count(mesh.cellPointIndices, sharedVertex)
        == static_cast<std::ranges::range_difference_t<decltype(mesh.cellPointIndices)>>(mesh.numberOfCells));

    hase::core::Point const origin = points[1u] * 0.2 + points[3u] * 0.1 + points[4u] * 0.3;
    hase::core::Point const direction = hase::kernels::forward::normalize(origin * -1.0);
    auto const result
        = hase::tests::traverseOneRay(device, backend[alpaka::object::exec], mesh, 0u, origin, direction);

    CHECK(hase::kernels::forward::hasMultipleTiedFaces(result.intersection.tiedFaceMask));
    CHECK(result.intersection.tiedFaceMask == ((1u << 1u) | (1u << 2u) | (1u << 3u)));
    CHECK(result.intersectionPoint.x == Catch::Approx(0.0).margin(1.0e-14));
    CHECK(result.intersectionPoint.y == Catch::Approx(0.0).margin(1.0e-14));
    CHECK(result.intersectionPoint.z == Catch::Approx(0.0).margin(1.0e-14));
    REQUIRE(result.transition.status == hase::kernels::forward::Tet4TransitionStatus::enteredCell);
    CHECK(result.transition.cell == 23u);
    CHECK(result.nextIntersection.localFace >= 0);
    CHECK(result.nextIntersection.length > 0.0);
}

TEMPLATE_LIST_TEST_CASE(
    "forward ray crosses a plane at perfectly perpendicular incidence",
    "[forward][traversal][edge-case]",
    hase::tests::TestBackends)
{
    auto const backend = TestType::makeDict();
    auto deviceSelector = alpaka::onHost::makeDeviceSelector(backend[alpaka::object::deviceSpec]);
    if(!deviceSelector.isAvailable())
    {
        SUCCEED("No device available for " << backend[alpaka::object::deviceSpec].getName());
        return;
    }
    auto device = deviceSelector.makeDevice(0);

    std::vector<hase::core::Point> const
        points{{-1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}};
    auto mesh = hase::tests::makeTraversalMesh(points, {{0u, 1u, 2u, 3u}, {4u, 1u, 3u, 2u}});

    hase::core::Point const direction{1.0, 0.0, 0.0};
    hase::core::Point const planeEdge0 = points[2u] - points[1u];
    hase::core::Point const planeEdge1 = points[3u] - points[1u];
    CHECK(hase::core::dot(direction, planeEdge0) == 0.0);
    CHECK(hase::core::dot(direction, planeEdge1) == 0.0);

    hase::core::Point const origin{-0.4, 0.2, 0.2};
    auto const result
        = hase::tests::traverseOneRay(device, backend[alpaka::object::exec], mesh, 0u, origin, direction);

    CHECK(result.intersection.localFace == 0);
    CHECK(result.intersection.tiedFaceMask == 1u);
    CHECK(result.intersection.length == Catch::Approx(0.4));
    CHECK(result.intersectionPoint.x == Catch::Approx(0.0).margin(1.0e-14));
    CHECK(result.intersectionPoint.y == Catch::Approx(0.2));
    CHECK(result.intersectionPoint.z == Catch::Approx(0.2));
    REQUIRE(result.transition.status == hase::kernels::forward::Tet4TransitionStatus::enteredCell);
    CHECK(result.transition.cell == 1u);
    CHECK(result.transition.forbiddenFace == 0);
    CHECK(result.nextIntersection.localFace >= 0);
    CHECK(result.nextIntersection.length == Catch::Approx(0.6));
}
