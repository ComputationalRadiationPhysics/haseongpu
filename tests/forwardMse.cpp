#include <alpaka/math.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <core/calcForwardPhiAse.hpp>
#include <kernels/forward/rayTransition.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <ranges>
#include <string>
#include <type_traits>
#include <vector>

namespace
{
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

    hase::core::DeviceMeshView traversalView(hase::core::HostMesh const& mesh)
    {
        hase::core::DeviceMeshView view{};
        view.points = mesh.points;
        view.cellPointIndices = mesh.cellPointIndices;
        view.cellFaces = mesh.cellFaces;
        view.barycentricFacePlanes = mesh.barycentricFacePlanes;
        view.cellNeighborCells = mesh.cellNeighborCells;
        view.cellNeighborLocalFaces = mesh.cellNeighborLocalFaces;
        view.numberOfCells = mesh.numberOfCells;
        view.numberOfFacesPerCell = mesh.numberOfFacesPerCell;
        view.numberOfCellVertices = mesh.numberOfCellVertices;
        view.numberOfMeshPoints = mesh.numberOfMeshPoints;
        return view;
    }
} // namespace

TEST_CASE("forward PhiASE RSE includes zero-score histories", "[forward][rse]")
{
    // Four globally launched histories with cell scores [1, 3, 0, 0].
    // The forward cell estimator scales the per-history mean by totalVolume / cellVolume.
    double const sum = 4.0;
    double const sumSquares = 10.0;
    unsigned const rayCount = 4u;
    double const totalVolume = 8.0;
    double const cellVolume = 4.0;

    double const expectedRelativeStandardError = std::sqrt((rayCount * sumSquares / (sum * sum) - 1.0) / rayCount);
    double const expectedStandardError = expectedRelativeStandardError * (sum * totalVolume / (rayCount * cellVolume));

    CHECK(
        hase::core::calcForwardRelativeStandardError(sum, sumSquares, rayCount)
        == Catch::Approx(expectedRelativeStandardError));
    CHECK(
        hase::core::calcForwardStandardError(sum, sumSquares, rayCount, totalVolume, cellVolume)
        == Catch::Approx(expectedStandardError));
}

TEST_CASE("forward PhiASE RSE handles invalid and zero-score estimates", "[forward][rse]")
{
    CHECK(hase::core::calcForwardRelativeStandardError(1.0, 1.0, 1u) == std::numeric_limits<double>::max());
    CHECK(alpaka::math::isnan(hase::core::calcForwardRelativeStandardError(0.0, 0.0, 2u)));
    CHECK(
        hase::core::calcForwardRelativeStandardError(std::numeric_limits<double>::infinity(), 1.0, 2u)
        == std::numeric_limits<double>::max());
    CHECK(hase::core::calcForwardStandardError(1.0, 1.0, 1u, 1.0, 1.0) == std::numeric_limits<double>::max());
    CHECK(hase::core::calcForwardStandardError(1.0, 1.0, 2u, 0.0, 1.0) == 0.0);
    CHECK(hase::core::calcForwardStandardError(0.0, 0.0, 2u, 1.0, 1.0) == 0.0);
    CHECK(hase::core::calcForwardStandardError(1.0, 1.0, 2u, 1.0, 0.0) == std::numeric_limits<double>::max());
    CHECK(
        hase::core::calcForwardStandardError(std::numeric_limits<double>::infinity(), 1.0, 2u, 1.0, 1.0)
        == std::numeric_limits<double>::max());
}

TEST_CASE("forward PhiASE beta-volume contribution uses double precision", "[forward][rse]")
{
    hase::core::BetaVolumeContribution contribution;
    auto const value = contribution(alpaka::Simd<double, 1u>{0.25}, alpaka::Simd<float, 1u>{0.5f});
    STATIC_REQUIRE(std::is_same_v<alpaka::trait::GetValueType_t<std::remove_cvref_t<decltype(value)>>, double>);
    CHECK(value[0] == Catch::Approx(0.125));
}

TEST_CASE("forward spectrum stratification balances discrete bins", "[forward][sampling]")
{
    constexpr unsigned spectrumSize = 7u;
    constexpr unsigned rayCount = 25u;
    std::array<unsigned, spectrumSize> visits{};
    for(unsigned ray = 0u; ray < rayCount; ++ray)
    {
        ++visits.at(hase::kernels::forward::stratifiedSpectrumIndex(spectrumSize, ray, rayCount, 3u));
    }

    auto const [minimum, maximum] = std::ranges::minmax_element(visits);
    CHECK(*maximum - *minimum <= 1u);
    CHECK(std::accumulate(visits.cbegin(), visits.cend(), 0u) == rayCount);
}

TEST_CASE("forward source stratification places one shifted point in each CDF interval", "[forward][sampling]")
{
    constexpr unsigned rayCount = 10u;
    constexpr double shift = 0.25;
    for(unsigned ray = 0u; ray < rayCount; ++ray)
    {
        double const target = hase::kernels::forward::stratifiedUnitInterval(ray, rayCount, shift);
        CHECK(target > static_cast<double>(ray) / rayCount);
        CHECK(target < static_cast<double>(ray + 1u) / rayCount);
    }
}

TEST_CASE("forward random histories are separated by ray, pass, and sampling domain", "[forward][sampling]")
{
    using hase::kernels::forward::rayHistoryId;
    using hase::kernels::forward::surfaceSamplingHistoryId;

    CHECK(rayHistoryId(0u, 7u) != rayHistoryId(0u, 8u));
    CHECK(rayHistoryId(0u, 7u) != rayHistoryId(1u, 7u));
    CHECK(rayHistoryId(1u, 7u) != surfaceSamplingHistoryId(1u));

    constexpr unsigned seed = 1234u;
    auto first = alpaka::rand::engine::Philox4x32x10{seed, rayHistoryId(3u, 11u)};
    auto repeated = alpaka::rand::engine::Philox4x32x10{seed, rayHistoryId(3u, 11u)};
    auto otherRay = alpaka::rand::engine::Philox4x32x10{seed, rayHistoryId(3u, 12u)};

    CHECK(first() == repeated());
    CHECK(first() != otherRay());
}

TEST_CASE("forward Tet4 face planes are barycentric", "[forward][traversal]")
{
    hase::core::HostMesh mesh;
    mesh.points = {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    mesh.numberOfCells = 1u;
    mesh.numberOfMeshPoints = 4u;
    mesh.cellPointIndices = {0u, 1u, 2u, 3u};
    mesh.cellFaces = {1, 2, 3, 0, 3, 2, 0, 1, 3, 0, 2, 1};
    mesh.precomputeBarycentricFacePlanes();

    auto const point = [](unsigned const vertex)
    {
        return std::array<hase::core::Point, 4u>{
            hase::core::Point{0.0, 0.0, 0.0},
            hase::core::Point{1.0, 0.0, 0.0},
            hase::core::Point{0.0, 1.0, 0.0},
            hase::core::Point{0.0, 0.0, 1.0}}
            .at(vertex);
    };
    auto const coordinate = [&mesh](unsigned const face, hase::core::Point const value)
    {
        unsigned const offset = face * hase::core::tet4BarycentricPlaneWidth;
        return mesh.barycentricFacePlanes[offset] * value.x + mesh.barycentricFacePlanes[offset + 1u] * value.y
               + mesh.barycentricFacePlanes[offset + 2u] * value.z + mesh.barycentricFacePlanes[offset + 3u];
    };

    for(unsigned face = 0u; face < hase::core::tet4FaceCount; ++face)
    {
        CHECK(coordinate(face, point(face)) == Catch::Approx(1.0));
        for(unsigned localVertex = 0u; localVertex < hase::core::tet4FaceWidth; ++localVertex)
        {
            CHECK(
                coordinate(
                    face,
                    point(static_cast<unsigned>(mesh.cellFaces[face * hase::core::tet4FaceWidth + localVertex])))
                == Catch::Approx(0.0));
        }
    }
    hase::core::Point const center{0.25, 0.25, 0.25};
    for(unsigned face = 0u; face < hase::core::tet4FaceCount; ++face)
        CHECK(coordinate(face, center) == Catch::Approx(0.25));

    CHECK(hase::kernels::forward::barycentricFaceIntersectionLength(0.3, -0.2, 2.0) == Catch::Approx(1.5));
    CHECK(hase::kernels::forward::barycentricFaceIntersectionLength(0.3, 0.2, 2.0) == 0.0);
    CHECK(hase::kernels::forward::barycentricFaceIntersectionLength(0.3, -0.2, 1.0) == 0.0);

    hase::core::DeviceMeshView view{};
    view.points = mesh.points;
    view.cellPointIndices = mesh.cellPointIndices;
    view.cellFaces = mesh.cellFaces;
    view.barycentricFacePlanes = mesh.barycentricFacePlanes;
    view.numberOfCells = 1u;
    view.numberOfFacesPerCell = hase::core::tet4FaceCount;
    view.numberOfCellVertices = hase::core::tet4VertexCount;
    view.numberOfMeshPoints = 4u;
    auto const intersection = hase::kernels::forward::nextFaceIntersection(
        view,
        0u,
        hase::core::Point{0.25, 0.25, 0.25},
        hase::core::Point{1.0, 0.0, 0.0},
        -1);
    CHECK(intersection.localFace == 0);
    CHECK(intersection.length == Catch::Approx(0.25));
    CHECK(intersection.tiedFaceMask == 1u);

    auto const twoFaceTie = hase::kernels::forward::nextFaceIntersection(
        view,
        0u,
        center,
        hase::kernels::forward::normalize(hase::core::Point{-1.0, -1.0, 0.0}),
        -1);
    CHECK(twoFaceTie.localFace == 1);
    CHECK(twoFaceTie.tiedFaceMask == ((1u << 1u) | (1u << 2u)));

    auto const threeFaceTie = hase::kernels::forward::nextFaceIntersection(
        view,
        0u,
        center,
        hase::kernels::forward::normalize(hase::core::Point{-1.0, -1.0, -1.0}),
        -1);
    CHECK(threeFaceTie.localFace == 1);
    CHECK(threeFaceTie.tiedFaceMask == ((1u << 1u) | (1u << 2u) | (1u << 3u)));

    hase::core::Point const reflectedOrigin = hase::kernels::forward::faceCentroid(view, 0u, 0u);
    auto const reflectedIntersection = hase::kernels::forward::nextFaceIntersection(
        view,
        0u,
        reflectedOrigin,
        hase::kernels::forward::normalize(hase::core::Point{-1.0, -2.0, -3.0}),
        0);
    CHECK(reflectedIntersection.localFace == 3);
    CHECK(reflectedIntersection.length > 0.0);
}

TEST_CASE("forward Tet4 intersection distances follow mesh scale", "[forward][traversal]")
{
    for(double const scale : {1.0e-9, 1.0, 1.0e9})
    {
        std::array<double, hase::core::tet4FaceCount * hase::core::tet4BarycentricPlaneWidth> planes{
            -1.0 / scale,
            -1.0 / scale,
            -1.0 / scale,
            1.0,
            1.0 / scale,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0 / scale,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0 / scale,
            0.0};
        hase::core::DeviceMeshView view{};
        view.barycentricFacePlanes = planes;
        view.numberOfCells = 1u;
        view.numberOfFacesPerCell = hase::core::tet4FaceCount;

        auto const intersection = hase::kernels::forward::nextFaceIntersection(
            view,
            0u,
            hase::core::Point{0.25 * scale, 0.25 * scale, 0.25 * scale},
            hase::core::Point{1.0, 0.0, 0.0},
            -1);
        CHECK(intersection.localFace == 0);
        CHECK(intersection.tiedFaceMask == 1u);
        CHECK(intersection.length == Catch::Approx(0.25 * scale));
    }
}

TEST_CASE("thin Tet4 mesh reproduces the old nudge lost-ray failure", "[forward][traversal]")
{
    constexpr double longCellHeight = 1.0e9;
    constexpr double thinCellHeight = 1.0e-7;
    auto const mesh = makeTraversalMesh(
        {
            {-longCellHeight, 0.0, 0.0},
            {0.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
            {thinCellHeight, 0.0, 0.0},
        },
        {{0u, 1u, 2u, 3u}, {4u, 1u, 2u, 3u}});
    auto const view = traversalView(mesh);
    hase::core::Point const direction{1.0, 0.0, 0.0};
    hase::core::Point const origin{-0.5 * longCellHeight, 0.125, 0.125};

    auto const longCellIntersection = hase::kernels::forward::nextFaceIntersection(view, 0u, origin, direction, -1);
    REQUIRE(longCellIntersection.localFace == 0);
    REQUIRE(longCellIntersection.tiedFaceMask == 1u);
    hase::core::Point const sharedFacePoint
        = hase::kernels::forward::advance(origin, direction, longCellIntersection.length);

    auto const transition = hase::kernels::forward::transitionAcrossIntersection(
        view,
        0u,
        longCellIntersection,
        sharedFacePoint,
        direction);
    REQUIRE(transition.status == hase::kernels::forward::Tet4TransitionStatus::enteredCell);
    REQUIRE(transition.cell == 1u);

    auto const exactThinCellIntersection = hase::kernels::forward::nextFaceIntersection(
        view,
        transition.cell,
        sharedFacePoint,
        direction,
        transition.forbiddenFace);
    REQUIRE(exactThinCellIntersection.localFace >= 0);
    CHECK(exactThinCellIntersection.length > 0.0);
    CHECK(exactThinCellIntersection.length < 2.0 * thinCellHeight);

    double const oldNudge = 64.0 * std::numeric_limits<double>::epsilon() * longCellIntersection.length;
    REQUIRE(oldNudge > exactThinCellIntersection.length);
    hase::core::Point const oldNudgedOrigin = hase::kernels::forward::advance(sharedFacePoint, direction, oldNudge);
    auto const lostRayIntersection
        = hase::kernels::forward::nextFaceIntersection(view, 1u, oldNudgedOrigin, direction, 0);
    // This was the walker's droppedRays branch: the nudge has skipped the thin cell completely.
    CHECK(lostRayIntersection.localFace < 0);
}

TEST_CASE("forward Tet4 recovery crosses a shared edge", "[forward][traversal]")
{
    auto const mesh = makeTraversalMesh(
        {
            {0.0, 0.0, -1.0},
            {0.0, 0.0, 1.0},
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {-1.0, 0.0, 0.0},
            {0.0, -1.0, 0.0},
        },
        {
            {0u, 1u, 2u, 3u},
            {0u, 1u, 3u, 4u},
            {0u, 1u, 4u, 5u},
            {0u, 1u, 5u, 2u},
        });
    auto const view = traversalView(mesh);
    hase::core::Point const origin{0.25, 0.25, 0.0};
    hase::core::Point const direction = hase::kernels::forward::normalize(hase::core::Point{-1.0, -1.0, 0.0});
    auto const intersection = hase::kernels::forward::nextFaceIntersection(view, 0u, origin, direction, -1);
    REQUIRE(hase::kernels::forward::hasMultipleTiedFaces(intersection.tiedFaceMask));
    hase::core::Point const edgePoint = hase::kernels::forward::advance(origin, direction, intersection.length);

    auto const transition
        = hase::kernels::forward::transitionAcrossIntersection(view, 0u, intersection, edgePoint, direction);
    CHECK(transition.status == hase::kernels::forward::Tet4TransitionStatus::enteredCell);
    CHECK(transition.cell == 2u);
}

TEST_CASE("forward Tet4 recovery crosses a shared vertex", "[forward][traversal]")
{
    auto const mesh = makeTraversalMesh(
        {
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {-1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, -1.0, 0.0},
            {0.0, 0.0, 1.0},
            {0.0, 0.0, -1.0},
        },
        {
            {0u, 1u, 3u, 5u},
            {0u, 2u, 3u, 5u},
            {0u, 2u, 4u, 5u},
            {0u, 2u, 4u, 6u},
        });
    auto const view = traversalView(mesh);
    hase::core::Point const origin{0.25, 0.25, 0.25};
    hase::core::Point const direction = hase::kernels::forward::normalize(hase::core::Point{-1.0, -1.0, -1.0});
    auto const intersection = hase::kernels::forward::nextFaceIntersection(view, 0u, origin, direction, -1);
    REQUIRE(hase::kernels::forward::hasMultipleTiedFaces(intersection.tiedFaceMask));
    hase::core::Point const vertex = hase::kernels::forward::advance(origin, direction, intersection.length);

    auto const transition
        = hase::kernels::forward::transitionAcrossIntersection(view, 0u, intersection, vertex, direction);
    CHECK(transition.status == hase::kernels::forward::Tet4TransitionStatus::enteredCell);
    CHECK(transition.cell == 3u);
}

TEST_CASE("forward Tet4 probe recovery selects an alternate neighbor", "[forward][traversal]")
{
    constexpr unsigned numberOfCells = 3u;
    constexpr unsigned faceCount = hase::core::tet4FaceCount;
    constexpr unsigned planeWidth = hase::core::tet4BarycentricPlaneWidth;
    std::array<double, numberOfCells * faceCount * planeWidth> planes{};
    auto const setNegativeProbeFace = [&planes](unsigned const cell, unsigned const face)
    { planes[(cell * faceCount + face) * planeWidth] = -1.0; };
    setNegativeProbeFace(0u, 0u);
    setNegativeProbeFace(0u, 1u);
    setNegativeProbeFace(1u, 0u);

    std::array<int, numberOfCells * faceCount> neighbors{};
    std::array<int, numberOfCells * faceCount> neighborFaces{};
    neighbors.fill(-1);
    neighborFaces.fill(-1);
    neighbors[0u * faceCount + 0u] = 1;
    neighborFaces[0u * faceCount + 0u] = 0;
    neighbors[0u * faceCount + 1u] = 2;
    neighborFaces[0u * faceCount + 1u] = 0;
    neighbors[1u * faceCount + 0u] = 0;
    neighborFaces[1u * faceCount + 0u] = 0;
    neighbors[2u * faceCount + 0u] = 0;
    neighborFaces[2u * faceCount + 0u] = 1;

    hase::core::DeviceMeshView view{};
    view.barycentricFacePlanes = planes;
    view.cellNeighborCells = neighbors;
    view.cellNeighborLocalFaces = neighborFaces;
    view.numberOfCells = numberOfCells;
    view.numberOfFacesPerCell = faceCount;

    hase::core::Point const hitPoint{0.0, 0.0, 0.0};
    auto const transition
        = hase::kernels::forward::recoverFaceTransition(view, 0u, 0, hitPoint, hase::core::Point{1.0, 0.0, 0.0});

    CHECK(transition.status == hase::kernels::forward::Tet4TransitionStatus::enteredCell);
    CHECK(transition.cell == 2u);
    CHECK(hitPoint.x == 0.0);
    CHECK(hitPoint.y == 0.0);
    CHECK(hitPoint.z == 0.0);
}

TEST_CASE("forward Tet4 recovery remains bounded on cyclic connectivity", "[forward][traversal]")
{
    std::array<double, hase::core::tet4FaceCount * hase::core::tet4BarycentricPlaneWidth> planes{};
    planes[0u] = -1.0;
    std::array<int, hase::core::tet4FaceCount> neighbors{0, -1, -1, -1};
    std::array<int, hase::core::tet4FaceCount> neighborFaces{1, -1, -1, -1};

    hase::core::DeviceMeshView view{};
    view.barycentricFacePlanes = planes;
    view.cellNeighborCells = neighbors;
    view.cellNeighborLocalFaces = neighborFaces;
    view.numberOfCells = 1u;
    view.numberOfFacesPerCell = hase::core::tet4FaceCount;

    auto const transition = hase::kernels::forward::recoverFaceTransition(
        view,
        0u,
        0,
        hase::core::Point{0.0, 0.0, 0.0},
        hase::core::Point{1.0, 0.0, 0.0});
    CHECK(transition.status == hase::kernels::forward::Tet4TransitionStatus::failed);
}

TEST_CASE("forward SRM environment controls are strict positive overrides", "[forward][srm]")
{
    auto const restore = [](char const* name, char const* value)
    {
        if(value == nullptr)
            unsetenv(name);
        else
            setenv(name, value, 1);
    };
    char const* oldMaxIterations = std::getenv("HASE_SRM_MAX_ITERATIONS");
    char const* oldDivergenceStreak = std::getenv("HASE_SRM_DIVERGENCE_STREAK");
    std::string const savedMaxIterations = oldMaxIterations == nullptr ? "" : oldMaxIterations;
    std::string const savedDivergenceStreak = oldDivergenceStreak == nullptr ? "" : oldDivergenceStreak;

    unsetenv("HASE_SRM_MAX_ITERATIONS");
    unsetenv("HASE_SRM_DIVERGENCE_STREAK");
    hase::core::ExperimentParameters experiment{};
    experiment.reflectionMaxIterations = 8u;
    auto const defaults = hase::core::resolveSrmControls(experiment);
    CHECK(defaults.maxIterations == 8u);
    CHECK(defaults.divergenceStreak == 3u);

    setenv("HASE_SRM_MAX_ITERATIONS", "11", 1);
    setenv("HASE_SRM_DIVERGENCE_STREAK", "4", 1);
    auto const overridden = hase::core::resolveSrmControls(experiment);
    CHECK(overridden.maxIterations == 11u);
    CHECK(overridden.divergenceStreak == 4u);

    setenv("HASE_SRM_DIVERGENCE_STREAK", "0", 1);
    CHECK_THROWS_AS(hase::core::resolveSrmControls(experiment), std::runtime_error);

    restore("HASE_SRM_MAX_ITERATIONS", oldMaxIterations == nullptr ? nullptr : savedMaxIterations.c_str());
    restore("HASE_SRM_DIVERGENCE_STREAK", oldDivergenceStreak == nullptr ? nullptr : savedDivergenceStreak.c_str());
}
