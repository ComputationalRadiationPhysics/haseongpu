/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <alpaka/alpaka.hpp>

#include <alpakaUtils/DevBundle.hpp>
#include <alpakaUtils/memory.hpp>
#include <alpakaUtils/utils.hpp>
#include <core/hostRoutineTiming.hpp>
#include <core/mesh.hpp>
#include <core/simulationRunControl.hpp>
#include <kernels/forward/barycentric.hpp>
#include <kernels/forward/policyRay.hpp>
#include <kernels/forward/rayTransition.hpp>
#include <kernels/forward/rayWalk.hpp>
#include <kernels/vertexAccumulation.hpp>
#include <random/random.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <iterator>
#include <limits>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace hase::kernels
{
    struct PumpBoundaryFace
    {
        unsigned cell = 0u;
        unsigned localFace = 0u;
        int domain = 0;
        std::array<hase::core::Point, 3u> vertices;
        hase::core::Point centroid;
        hase::core::Point normal;
        double area = 0.0;
    };

    struct GeneralPumpRayState
        : hase::kernels::forward::ray::TraversalState
        , hase::kernels::forward::ray::SrmPositionStorage<
              typename std::remove_cvref_t<ALPAKA_TYPEOF(hase::kernels::forward::ray::pumpSrmPolicy)>::PositionPolicy>
    {
        double power = 0.0;
        double wavelength = 0.0;
        double sigmaAbsorption = 0.0;
        double sigmaEmission = 0.0;
        unsigned relayIndex = 0u;
    };

    static_assert(std::derived_from<GeneralPumpRayState, hase::kernels::forward::ray::BarycentricSrmPositionStorage>);

    [[nodiscard]] inline bool containsDomain(std::vector<int> const& domains, int const domain)
    {
        return std::find(domains.begin(), domains.end(), domain) != domains.end();
    }

    [[nodiscard]] inline hase::core::Point hostPoint(hase::core::HostMesh const& mesh, unsigned const point)
    {
        return {
            mesh.points[point],
            mesh.points[point + mesh.numberOfMeshPoints],
            mesh.points[point + 2u * mesh.numberOfMeshPoints]};
    }

    [[nodiscard]] inline std::vector<PumpBoundaryFace> pumpBoundaryFaces(
        hase::core::HostMesh const& mesh,
        std::vector<int> const& domains)
    {
        std::vector<PumpBoundaryFace> result;
        for(unsigned cell = 0u; cell < mesh.numberOfCells; ++cell)
        {
            for(unsigned face = 0u; face < mesh.numberOfFacesPerCell; ++face)
            {
                unsigned const faceIndex = cell * mesh.numberOfFacesPerCell + face;
                int const domain = mesh.cellFaceBoundaries[faceIndex];
                if(mesh.cellNeighborCells[faceIndex] >= 0 || !containsDomain(domains, domain))
                    continue;
                PumpBoundaryFace info;
                info.cell = cell;
                info.localFace = face;
                info.domain = domain;
                for(unsigned vertex = 0u; vertex < 3u; ++vertex)
                {
                    int const point = mesh.cellFaces[faceIndex * 3u + vertex];
                    if(point < 0)
                        throw std::runtime_error("pump boundary face contains an invalid point");
                    info.vertices[vertex] = hostPoint(mesh, static_cast<unsigned>(point));
                }
                info.centroid = (info.vertices[0] + info.vertices[1] + info.vertices[2]) * (1.0 / 3.0);
                auto normal
                    = hase::core::cross(info.vertices[1] - info.vertices[0], info.vertices[2] - info.vertices[0]);
                double const twiceArea = normal.euclidLength();
                if(twiceArea <= 0.0)
                    continue;
                info.area = 0.5 * twiceArea;
                info.normal = normal * (1.0 / twiceArea);
                hase::core::Point const center{
                    mesh.cellCenters[cell],
                    mesh.cellCenters[cell + mesh.numberOfCells],
                    mesh.cellCenters[cell + 2u * mesh.numberOfCells]};
                if(hase::core::dot(info.normal, center - info.centroid) > 0.0)
                    info.normal = info.normal * -1.0;
                result.push_back(info);
            }
        }
        return result;
    }

    [[nodiscard]] inline hase::core::Point hostNormalize(hase::core::Point const value)
    {
        double const length = value.euclidLength();
        if(length <= 0.0)
            return {0.0, 0.0, 0.0};
        return value * (1.0 / length);
    }

    [[nodiscard]] inline hase::core::Point perpendicular(hase::core::Point const normal)
    {
        hase::core::Point reference
            = std::abs(normal.x) < 0.9 ? hase::core::Point{1.0, 0.0, 0.0} : hase::core::Point{0.0, 1.0, 0.0};
        return hostNormalize(hase::core::cross(normal, reference));
    }

    [[nodiscard]] inline double pumpProfileWeight(
        hase::core::PumpProfileParameters const& profile,
        hase::core::Point const point)
    {
        if(profile.kind == 0u)
            return 1.0;
        hase::core::Point const relative
            = point - hase::core::Point{profile.center[0], profile.center[1], profile.center[2]};
        double const u
            = hase::core::dot(relative, hase::core::Point{profile.axisU[0], profile.axisU[1], profile.axisU[2]})
              / profile.radiusU;
        double const v
            = hase::core::dot(relative, hase::core::Point{profile.axisV[0], profile.axisV[1], profile.axisV[2]})
              / profile.radiusV;
        return std::exp(-std::pow(std::sqrt(u * u + v * v), profile.exponent));
    }

    [[nodiscard]] inline double pumpEntryWeight(
        PumpBoundaryFace const& face,
        hase::core::PumpProfileParameters const& profile)
    {
        // Seven-point Dunavant quadrature.  The CDF must represent the
        // aperture's spatial entry distribution, not just its triangle areas:
        // after a region is selected, rejection sampling supplies the matching
        // conditional distribution inside that region.
        constexpr std::array<std::array<double, 3u>, 7u> barycentric{{
            {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
            {0.059715871789770, 0.470142064105115, 0.470142064105115},
            {0.470142064105115, 0.059715871789770, 0.470142064105115},
            {0.470142064105115, 0.470142064105115, 0.059715871789770},
            {0.797426985353087, 0.101286507323456, 0.101286507323456},
            {0.101286507323456, 0.797426985353087, 0.101286507323456},
            {0.101286507323456, 0.101286507323456, 0.797426985353087},
        }};
        constexpr std::array<double, 7u> weights{
            0.225,
            0.132394152788506,
            0.132394152788506,
            0.132394152788506,
            0.125939180544827,
            0.125939180544827,
            0.125939180544827};

        double integral = 0.0;
        for(std::size_t sample = 0u; sample < barycentric.size(); ++sample)
        {
            auto const& coordinate = barycentric[sample];
            hase::core::Point const point = face.vertices[0] * coordinate[0] + face.vertices[1] * coordinate[1]
                                            + face.vertices[2] * coordinate[2];
            integral += weights[sample] * pumpProfileWeight(profile, point);
        }
        return face.area * integral;
    }

    inline void appendPumpEntryRegions(
        PumpBoundaryFace const& face,
        unsigned const remainingSubdivisions,
        std::vector<PumpBoundaryFace>& regions)
    {
        if(remainingSubdivisions == 0u)
        {
            regions.push_back(face);
            return;
        }

        auto const midpoint = [](hase::core::Point const a, hase::core::Point const b) { return (a + b) * 0.5; };
        hase::core::Point const midpoint01 = midpoint(face.vertices[0], face.vertices[1]);
        hase::core::Point const midpoint12 = midpoint(face.vertices[1], face.vertices[2]);
        hase::core::Point const midpoint20 = midpoint(face.vertices[2], face.vertices[0]);
        std::array<std::array<hase::core::Point, 3u>, 4u> const vertices{{
            {face.vertices[0], midpoint01, midpoint20},
            {midpoint01, face.vertices[1], midpoint12},
            {midpoint20, midpoint12, face.vertices[2]},
            {midpoint01, midpoint12, midpoint20},
        }};
        for(auto const& regionVertices : vertices)
        {
            PumpBoundaryFace region = face;
            region.vertices = regionVertices;
            region.centroid = (regionVertices[0] + regionVertices[1] + regionVertices[2]) * (1.0 / 3.0);
            region.area = face.area * 0.25;
            appendPumpEntryRegions(region, remainingSubdivisions - 1u, regions);
        }
    }

    [[nodiscard]] inline std::vector<PumpBoundaryFace> pumpEntryRegions(std::vector<PumpBoundaryFace> const& faces)
    {
        // Spatial regions make the systematic CDF cover the continuous entry
        // aperture.  A face-only CDF leaves all within-face variation random,
        // which is especially noisy for a nonuniform beam profile.
        constexpr unsigned subdivisionDepth = 2u;
        constexpr std::size_t regionsPerFace = 1u << (2u * subdivisionDepth);
        std::vector<PumpBoundaryFace> regions;
        regions.reserve(faces.size() * regionsPerFace);
        for(auto const& face : faces)
            appendPumpEntryRegions(face, subdivisionDepth, regions);
        return regions;
    }

    [[nodiscard]] inline std::size_t pumpEntryRegionForTarget(std::vector<double> const& entryCdf, double const target)
    {
        if(entryCdf.empty())
            return 0u;
        auto const region = std::upper_bound(entryCdf.cbegin(), entryCdf.cend(), target);
        return region == entryCdf.cend() ? entryCdf.size() - 1u
                                         : static_cast<std::size_t>(std::distance(entryCdf.cbegin(), region));
    }

    template<typename T_Rng>
    [[nodiscard]] inline hase::core::Point sampleTriangle(PumpBoundaryFace const& face, T_Rng& rng)
    {
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        double u = uniform(rng);
        double v = uniform(rng);
        if(u + v > 1.0)
        {
            u = 1.0 - u;
            v = 1.0 - v;
        }
        return face.vertices[0] + (face.vertices[1] - face.vertices[0]) * u
               + (face.vertices[2] - face.vertices[0]) * v;
    }

    [[nodiscard]] inline std::vector<GeneralPumpRayState> samplePumpSource(
        hase::core::HostMesh const& mesh,
        hase::core::PumpSourceParameters const& source,
        unsigned const globalRayCount,
        std::uint32_t const seed,
        unsigned const firstRay = 0u,
        unsigned const localRayCount = std::numeric_limits<unsigned>::max())
    {
        unsigned const selectedRayCount = std::min(localRayCount, globalRayCount - std::min(firstRay, globalRayCount));
        unsigned const selectedRayEnd = std::min(globalRayCount, firstRay + selectedRayCount);
        auto const faces = pumpBoundaryFaces(mesh, source.surfaces);
        if(faces.empty())
            throw std::runtime_error("pump source selected no exterior boundary faces");
        if(selectedRayCount == 0u)
            return {};
        auto const entryRegions = pumpEntryRegions(faces);
        std::vector<double> entryCdf;
        entryCdf.reserve(entryRegions.size());
        double totalEntryWeight = 0.0;
        for(auto const& region : entryRegions)
        {
            totalEntryWeight += pumpEntryWeight(region, source.profile);
            entryCdf.push_back(totalEntryWeight);
        }
        if(!(totalEntryWeight > 0.0) || !std::isfinite(totalEntryWeight))
            throw std::runtime_error("pump spatial profile has no finite weight on the selected exterior faces");
        std::discrete_distribution<std::size_t> spectrumDistribution(
            source.spectralWeights.begin(),
            source.spectralWeights.end());
        std::discrete_distribution<std::size_t> angularDistribution(
            source.angularWeights.begin(),
            source.angularWeights.end());
        std::mt19937_64 rng(seed);
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        double const entryStratificationOffset = hase::random::stratifiedUnitOffset(seed);

        std::vector<GeneralPumpRayState> rays;
        rays.reserve(selectedRayCount);
        for(unsigned ray = 0u; ray < selectedRayEnd; ++ray)
        {
            double const entryTarget = (static_cast<double>(ray) + entryStratificationOffset)
                                       / static_cast<double>(globalRayCount) * totalEntryWeight;
            PumpBoundaryFace const* face = &entryRegions[pumpEntryRegionForTarget(entryCdf, entryTarget)];
            hase::core::Point origin;
            bool accepted = false;
            for(unsigned attempt = 0u; attempt < 100000u; ++attempt)
            {
                origin = sampleTriangle(*face, rng);
                if(uniform(rng) <= pumpProfileWeight(source.profile, origin))
                {
                    accepted = true;
                    break;
                }
            }
            if(!accepted)
                throw std::runtime_error("pump spatial profile rejection sampling did not converge");

            std::size_t const angular = angularDistribution(rng);
            double const theta = source.polarAngles[angular];
            double const phi = source.azimuthalAngles[angular];
            hase::core::Point const inward = face->normal * -1.0;
            hase::core::Point const u = perpendicular(inward);
            hase::core::Point const v = hase::core::cross(inward, u);
            hase::core::Point const direction = hostNormalize(
                inward * std::cos(theta) + u * (std::sin(theta) * std::cos(phi))
                + v * (std::sin(theta) * std::sin(phi)));
            std::size_t const spectrum = spectrumDistribution(rng);

            if(ray < firstRay)
                continue;

            GeneralPumpRayState rayState;
            rayState.position = origin;
            rayState.direction = direction;
            rayState.power = source.totalPower / static_cast<double>(globalRayCount);
            rayState.wavelength = source.wavelengths[spectrum];
            rayState.sigmaAbsorption = source.sigmaAbsorption[spectrum];
            rayState.sigmaEmission = source.sigmaEmission[spectrum];
            rayState.cell = face->cell;
            rayState.forbiddenFace = static_cast<int>(face->localFace);
            rays.push_back(rayState);
        }
        return rays;
    }

    struct StorePumpSrmBoundary
        : hase::kernels::forward::ray::BoundaryPolicySrm<hase::kernels::forward::ray::srmPosition::Barycentric>
    {
        ALPAKA_FN_ACC hase::kernels::forward::ray::BoundaryResult operator()(
            auto const&,
            hase::core::DeviceMeshView const& mesh,
            auto& ray,
            unsigned const cell,
            unsigned const localFace)
        {
            namespace policyRay = hase::kernels::forward::ray;
            policyRay::captureSrmPosition(this->positionPolicy, mesh, cell, localFace, ray.position, ray);
            return policyRay::BoundaryResult::stop;
        }
    };

    struct StorePumpSrmBoundaryFactory
    {
        ALPAKA_FN_ACC auto operator()(unsigned) const
        {
            return StorePumpSrmBoundary{};
        }
    };

    struct TraceGeneralPump
    {
        double planckConstant = 6.62607015e-34;
        double speedOfLight = 299792458.0;

        template<typename T_BetaVolumeView, typename T_VertexPumpIntegralView>
        struct CellPolicy : hase::kernels::forward::ray::behaviourDimension::Cell
        {
            double planckConstant;
            double speedOfLight;
            T_BetaVolumeView betaVolume;
            T_VertexPumpIntegralView vertexPumpIntegral;

            ALPAKA_FN_HOST_ACC constexpr CellPolicy(
                double const planckConstantValue,
                double const speedOfLightValue,
                T_BetaVolumeView betaVolumeValue,
                T_VertexPumpIntegralView vertexPumpIntegralValue)
                : planckConstant{planckConstantValue}
                , speedOfLight{speedOfLightValue}
                , betaVolume{betaVolumeValue}
                , vertexPumpIntegral{vertexPumpIntegralValue}
            {
            }

            ALPAKA_FN_ACC bool operator()(
                auto const& acc,
                hase::core::DeviceMeshView const& mesh,
                auto& ray,
                unsigned const tet,
                hase::kernels::forward::Tet4FaceIntersection const intersection)
            {
                bool const gainCell = mesh.getCellType(tet) != mesh.claddingNumber;
                double const gain
                    = gainCell
                          ? static_cast<double>(mesh.nTot)
                                * (betaVolume[tet] * (ray.sigmaAbsorption + ray.sigmaEmission) - ray.sigmaAbsorption)
                          : -mesh.claddingAbsorption;
                double const exponent = gain * intersection.length;
                if(!alpaka::math::isfinite(exponent) || exponent > 700.0)
                {
                    ray.power = 0.0;
                    return false;
                }
                double const nextPower = ray.power * alpaka::math::exp(exponent);
                if(gainCell && mesh.nTot > 0.0f)
                {
                    double const integral = (ray.power - nextPower) * ray.wavelength
                                            / (planckConstant * speedOfLight * static_cast<double>(mesh.nTot));
                    // Clamping and renormalizing protects positivity and exact integral
                    // conservation against round-off at faces.
                    auto const weights = hase::kernels::forward::segmentMidpointBarycentricVertexWeights(
                        mesh,
                        tet,
                        ray.position,
                        ray.direction,
                        intersection.length);
                    for(unsigned localVertex = 0u; localVertex < hase::core::tet4VertexCount; ++localVertex)
                    {
                        unsigned const point = mesh.cellPointIndices[tet * mesh.numberOfCellVertices + localVertex];
                        alpaka::onAcc::atomicAdd(acc, &vertexPumpIntegral[point], integral * weights[localVertex]);
                    }
                }
                ray.power = nextPower;
                return ray.power != 0.0;
            }
        };

        template<
            typename T_Acc,
            typename T_BetaVolumeView,
            typename T_OriginXView,
            typename T_OriginYView,
            typename T_OriginZView,
            typename T_DirectionXView,
            typename T_DirectionYView,
            typename T_DirectionZView,
            typename T_PowerView,
            typename T_WavelengthView,
            typename T_SigmaAbsorptionView,
            typename T_SigmaEmissionView,
            typename T_CellView,
            typename T_ForbiddenFaceView,
            typename T_BoundaryPolicyFactory,
            typename T_VertexPumpIntegralView>
        ALPAKA_FN_ACC void operator()(
            T_Acc const& acc,
            hase::core::DeviceMeshView const mesh,
            T_BetaVolumeView betaVolume,
            T_OriginXView originX,
            T_OriginYView originY,
            T_OriginZView originZ,
            T_DirectionXView directionX,
            T_DirectionYView directionY,
            T_DirectionZView directionZ,
            T_PowerView power,
            T_WavelengthView wavelength,
            T_SigmaAbsorptionView sigmaAbsorption,
            T_SigmaEmissionView sigmaEmission,
            T_CellView cell,
            T_ForbiddenFaceView forbiddenFace,
            T_BoundaryPolicyFactory boundaryPolicyFactory,
            T_VertexPumpIntegralView vertexPumpIntegral,
            unsigned const rayCount) const
        {
            for(auto [rayIndex] :
                alpaka::onAcc::makeIdxMap(acc, alpaka::onAcc::worker::threadsInGrid, alpaka::IdxRange{rayCount}))
            {
                namespace ray = hase::kernels::forward::ray;
                GeneralPumpRayState rayState;
                rayState.position = {originX[rayIndex], originY[rayIndex], originZ[rayIndex]};
                rayState.direction = {directionX[rayIndex], directionY[rayIndex], directionZ[rayIndex]};
                rayState.cell = cell[rayIndex];
                rayState.forbiddenFace = forbiddenFace[rayIndex];
                rayState.power = power[rayIndex];
                rayState.wavelength = wavelength[rayIndex];
                rayState.sigmaAbsorption = sigmaAbsorption[rayIndex];
                rayState.sigmaEmission = sigmaEmission[rayIndex];
                auto boundaryPolicy = boundaryPolicyFactory(rayIndex);
                static_cast<void>(ray::walk(
                    acc,
                    mesh,
                    rayState,
                    ray::RayWalkBehaviour{
                        CellPolicy<T_BetaVolumeView, T_VertexPumpIntegralView>{
                            planckConstant,
                            speedOfLight,
                            betaVolume,
                            vertexPumpIntegral},
                        boundaryPolicy}));
            }
        }
    };

    struct RelayFrame
    {
        hase::core::Point origin, u, v, normal;
        std::vector<PumpBoundaryFace> faces;
    };

    struct PumpRelayDeviceDescriptor
    {
        hase::core::Point exitOrigin, exitU, exitV, exitNormal;
        hase::core::Point entryOrigin, entryU, entryV, entryNormal;
        double cosine = 1.0;
        double sine = 0.0;
        double offsetU = 0.0;
        double offsetV = 0.0;
        double tiltU = 0.0;
        double tiltV = 0.0;
        double magnification = 1.0;
        //! Relay throughput applied before inward reinjection at the entry surface.
        double transmission = 1.0;
        int flipU = 1;
        int flipV = 1;
        unsigned entryFaceBegin = 0u;
        unsigned entryFaceEnd = 0u;
    };

    namespace pumpBoundaryPolicy
    {
        struct Relay
        {
        };

        struct SrmBarycentric
        {
        };
    } // namespace pumpBoundaryPolicy

    inline constexpr pumpBoundaryPolicy::Relay pumpRelayPolicy{};
    inline constexpr pumpBoundaryPolicy::SrmBarycentric pumpSrmBarycentricPolicy{};

    template<typename T_DescriptorView, typename T_UnsignedView, typename T_DoubleView>
    struct DevicePumpRelayBoundary
        : hase::kernels::forward::ray::BoundaryPolicySrm<hase::kernels::forward::ray::srmPosition::Barycentric>
    {
        T_DescriptorView descriptors;
        T_UnsignedView exitMask;
        T_UnsignedView entryFaceIds;
        T_UnsignedView cacheState;
        T_UnsignedView cacheTargetFace;
        T_DoubleView cacheBarycentric0;
        T_DoubleView cacheBarycentric1;
        T_DoubleView cacheBarycentric2;
        unsigned faceCount;
        unsigned relayCount;
        unsigned rayCount;
        unsigned rayIndex;

        ALPAKA_FN_HOST_ACC constexpr DevicePumpRelayBoundary(
            T_DescriptorView descriptorsValue,
            T_UnsignedView exitMaskValue,
            T_UnsignedView entryFaceIdsValue,
            T_UnsignedView cacheStateValue,
            T_UnsignedView cacheTargetFaceValue,
            T_DoubleView cacheBarycentric0Value,
            T_DoubleView cacheBarycentric1Value,
            T_DoubleView cacheBarycentric2Value,
            unsigned const faceCountValue,
            unsigned const relayCountValue,
            unsigned const rayCountValue,
            unsigned const rayIndexValue)
            : descriptors{descriptorsValue}
            , exitMask{exitMaskValue}
            , entryFaceIds{entryFaceIdsValue}
            , cacheState{cacheStateValue}
            , cacheTargetFace{cacheTargetFaceValue}
            , cacheBarycentric0{cacheBarycentric0Value}
            , cacheBarycentric1{cacheBarycentric1Value}
            , cacheBarycentric2{cacheBarycentric2Value}
            , faceCount{faceCountValue}
            , relayCount{relayCountValue}
            , rayCount{rayCountValue}
            , rayIndex{rayIndexValue}
        {
        }

        ALPAKA_FN_ACC hase::kernels::forward::ray::BoundaryResult operator()(
            auto const&,
            hase::core::DeviceMeshView const& mesh,
            auto& ray,
            unsigned const cell,
            unsigned const localFace)
        {
            namespace policyRay = hase::kernels::forward::ray;
            if(ray.relayIndex >= relayCount)
            {
                policyRay::captureSrmPosition(this->positionPolicy, mesh, cell, localFace, ray.position, ray);
                return policyRay::BoundaryResult::stop;
            }

            unsigned const relayIndex = ray.relayIndex;
            unsigned const faceId = cell * mesh.numberOfFacesPerCell + localFace;
            if(exitMask[relayIndex * faceCount + faceId] == 0u)
            {
                ray.power = 0.0;
                return policyRay::BoundaryResult::stop;
            }

            auto const& descriptor = descriptors[relayIndex];
            policyRay::captureSrmPosition(this->positionPolicy, mesh, cell, localFace, ray.position, ray);
            hase::core::Point const exitPosition
                = policyRay::restoreSrmPosition(this->positionPolicy, mesh, cell, localFace, ray);
            hase::core::Point const relative = exitPosition - descriptor.exitOrigin;
            double u = hase::core::dot(relative, descriptor.exitU) * static_cast<double>(descriptor.flipU);
            double v = hase::core::dot(relative, descriptor.exitV) * static_cast<double>(descriptor.flipV);
            u *= descriptor.magnification;
            v *= descriptor.magnification;
            double const mappedU = descriptor.cosine * u - descriptor.sine * v + descriptor.offsetU;
            double const mappedV = descriptor.sine * u + descriptor.cosine * v + descriptor.offsetV;
            hase::core::Point mappedPosition
                = descriptor.entryOrigin + descriptor.entryU * mappedU + descriptor.entryV * mappedV;

            unsigned const cacheIndex = relayIndex * rayCount + rayIndex;
            unsigned targetFaceId = cacheTargetFace[cacheIndex];
            policyRay::TriangleBarycentric targetCoordinates{
                cacheBarycentric0[cacheIndex],
                cacheBarycentric1[cacheIndex],
                cacheBarycentric2[cacheIndex]};
            unsigned const cachedState = cacheState[cacheIndex];
            if(cachedState == 0u)
            {
                targetFaceId = faceCount;
                for(unsigned entry = descriptor.entryFaceBegin; entry < descriptor.entryFaceEnd; ++entry)
                {
                    unsigned const candidateFaceId = entryFaceIds[entry];
                    unsigned const candidateCell = candidateFaceId / mesh.numberOfFacesPerCell;
                    unsigned const candidateLocalFace = candidateFaceId % mesh.numberOfFacesPerCell;
                    auto const coordinates = policyRay::triangleBarycentricCoordinates(
                        mesh,
                        candidateCell,
                        candidateLocalFace,
                        mappedPosition);
                    constexpr double tolerance = 1.0e-10;
                    if(coordinates[0u] >= -tolerance && coordinates[1u] >= -tolerance && coordinates[2u] >= -tolerance)
                    {
                        targetFaceId = candidateFaceId;
                        targetCoordinates = coordinates;
                        break;
                    }
                }
                cacheTargetFace[cacheIndex] = targetFaceId;
                cacheBarycentric0[cacheIndex] = targetCoordinates[0u];
                cacheBarycentric1[cacheIndex] = targetCoordinates[1u];
                cacheBarycentric2[cacheIndex] = targetCoordinates[2u];
                cacheState[cacheIndex] = targetFaceId < faceCount ? 1u : 2u;
            }
            else if(cachedState == 2u)
            {
                ray.power = 0.0;
                return policyRay::BoundaryResult::stop;
            }

            if(targetFaceId >= faceCount)
            {
                ray.power = 0.0;
                return policyRay::BoundaryResult::stop;
            }

            unsigned const targetCell = targetFaceId / mesh.numberOfFacesPerCell;
            unsigned const targetLocalFace = targetFaceId % mesh.numberOfFacesPerCell;
            mappedPosition
                = policyRay::positionFromTriangleBarycentric(mesh, targetCell, targetLocalFace, targetCoordinates);
            hase::core::Point const oldDirection = ray.direction;
            double const du = hase::core::dot(oldDirection, descriptor.exitU) * static_cast<double>(descriptor.flipU);
            double const dv = hase::core::dot(oldDirection, descriptor.exitV) * static_cast<double>(descriptor.flipV);
            double const mappedDu = descriptor.cosine * du - descriptor.sine * dv + descriptor.tiltU;
            double const mappedDv = descriptor.sine * du + descriptor.cosine * dv + descriptor.tiltV;
            double const normalMagnitude = alpaka::math::abs(hase::core::dot(oldDirection, descriptor.exitNormal));
            ray.direction = hase::kernels::forward::normalize(
                descriptor.entryU * mappedDu + descriptor.entryV * mappedDv
                - descriptor.entryNormal * normalMagnitude);
            ray.position = mappedPosition;
            ray.cell = targetCell;
            ray.forbiddenFace = static_cast<int>(targetLocalFace);
            ray.boundaryBarycentric = targetCoordinates;
            ray.power *= descriptor.transmission;
            ++ray.relayIndex;
            return ray.power == 0.0 ? policyRay::BoundaryResult::stop : policyRay::BoundaryResult::continueTraversal;
        }
    };

    template<typename T_DescriptorView, typename T_UnsignedView, typename T_DoubleView>
    struct DevicePumpRelayBoundaryFactory
    {
        T_DescriptorView descriptors;
        T_UnsignedView exitMask;
        T_UnsignedView entryFaceIds;
        T_UnsignedView cacheState;
        T_UnsignedView cacheTargetFace;
        T_DoubleView cacheBarycentric0;
        T_DoubleView cacheBarycentric1;
        T_DoubleView cacheBarycentric2;
        unsigned faceCount;
        unsigned relayCount;
        unsigned rayCount;

        ALPAKA_FN_ACC auto operator()(unsigned const rayIndex) const
        {
            return DevicePumpRelayBoundary{
                descriptors,
                exitMask,
                entryFaceIds,
                cacheState,
                cacheTargetFace,
                cacheBarycentric0,
                cacheBarycentric1,
                cacheBarycentric2,
                faceCount,
                relayCount,
                rayCount,
                rayIndex};
        }
    };

    template<typename T_DescriptorView, typename T_UnsignedView, typename T_DoubleView>
    struct PumpBarycentricSrmBoundary
        : hase::kernels::forward::ray::BoundaryPolicySrm<hase::kernels::forward::ray::srmPosition::Barycentric>
    {
        T_DescriptorView descriptors;
        T_UnsignedView surfaceMask;
        T_UnsignedView cacheState;
        T_UnsignedView cacheTargetFace;
        T_DoubleView cacheBarycentric0;
        T_DoubleView cacheBarycentric1;
        T_DoubleView cacheBarycentric2;
        unsigned faceCount;
        unsigned reflectionCount;
        unsigned rayCount;
        unsigned rayIndex;

        ALPAKA_FN_HOST_ACC constexpr PumpBarycentricSrmBoundary(
            T_DescriptorView descriptorsValue,
            T_UnsignedView surfaceMaskValue,
            T_UnsignedView cacheStateValue,
            T_UnsignedView cacheTargetFaceValue,
            T_DoubleView cacheBarycentric0Value,
            T_DoubleView cacheBarycentric1Value,
            T_DoubleView cacheBarycentric2Value,
            unsigned const faceCountValue,
            unsigned const reflectionCountValue,
            unsigned const rayCountValue,
            unsigned const rayIndexValue)
            : descriptors{descriptorsValue}
            , surfaceMask{surfaceMaskValue}
            , cacheState{cacheStateValue}
            , cacheTargetFace{cacheTargetFaceValue}
            , cacheBarycentric0{cacheBarycentric0Value}
            , cacheBarycentric1{cacheBarycentric1Value}
            , cacheBarycentric2{cacheBarycentric2Value}
            , faceCount{faceCountValue}
            , reflectionCount{reflectionCountValue}
            , rayCount{rayCountValue}
            , rayIndex{rayIndexValue}
        {
        }

        ALPAKA_FN_ACC hase::kernels::forward::ray::BoundaryResult operator()(
            auto const&,
            hase::core::DeviceMeshView const& mesh,
            auto& ray,
            unsigned const cell,
            unsigned const localFace)
        {
            namespace policyRay = hase::kernels::forward::ray;
            if(ray.relayIndex >= reflectionCount)
                return policyRay::BoundaryResult::stop;

            unsigned const reflectionIndex = ray.relayIndex;
            unsigned const faceId = cell * mesh.numberOfFacesPerCell + localFace;
            if(surfaceMask[reflectionIndex * faceCount + faceId] == 0u)
            {
                ray.power = 0.0;
                return policyRay::BoundaryResult::stop;
            }

            unsigned const cacheIndex = reflectionIndex * rayCount + rayIndex;
            if(cacheState[cacheIndex] == 0u)
            {
                policyRay::captureSrmPosition(this->positionPolicy, mesh, cell, localFace, ray.position, ray);
                cacheTargetFace[cacheIndex] = faceId;
                cacheBarycentric0[cacheIndex] = ray.boundaryBarycentric[0u];
                cacheBarycentric1[cacheIndex] = ray.boundaryBarycentric[1u];
                cacheBarycentric2[cacheIndex] = ray.boundaryBarycentric[2u];
                cacheState[cacheIndex] = 1u;
            }
            else
            {
                unsigned const cachedFace = cacheTargetFace[cacheIndex];
                if(cachedFace != faceId)
                {
                    ray.power = 0.0;
                    return policyRay::BoundaryResult::stop;
                }
                ray.boundaryBarycentric
                    = {cacheBarycentric0[cacheIndex], cacheBarycentric1[cacheIndex], cacheBarycentric2[cacheIndex]};
                ray.position = policyRay::restoreSrmPosition(this->positionPolicy, mesh, cell, localFace, ray);
            }

            ray.direction = hase::kernels::forward::reflectedDirection(
                ray.direction,
                hase::kernels::forward::outwardFaceNormal(mesh, cell, localFace));
            ray.cell = cell;
            ray.forbiddenFace = static_cast<int>(localFace);
            ray.power *= descriptors[reflectionIndex].transmission;
            ++ray.relayIndex;
            return ray.power == 0.0 ? policyRay::BoundaryResult::stop : policyRay::BoundaryResult::continueTraversal;
        }
    };

    template<typename T_DescriptorView, typename T_UnsignedView, typename T_DoubleView>
    struct PumpBarycentricSrmBoundaryFactory
    {
        T_DescriptorView descriptors;
        T_UnsignedView surfaceMask;
        T_UnsignedView cacheState;
        T_UnsignedView cacheTargetFace;
        T_DoubleView cacheBarycentric0;
        T_DoubleView cacheBarycentric1;
        T_DoubleView cacheBarycentric2;
        unsigned faceCount;
        unsigned reflectionCount;
        unsigned rayCount;

        ALPAKA_FN_ACC auto operator()(unsigned const rayIndex) const
        {
            return PumpBarycentricSrmBoundary{
                descriptors,
                surfaceMask,
                cacheState,
                cacheTargetFace,
                cacheBarycentric0,
                cacheBarycentric1,
                cacheBarycentric2,
                faceCount,
                reflectionCount,
                rayCount,
                rayIndex};
        }
    };

    [[nodiscard]] inline RelayFrame makeRelayFrame(hase::core::HostMesh const& mesh, std::vector<int> const& domains)
    {
        RelayFrame frame{};
        frame.faces = pumpBoundaryFaces(mesh, domains);
        if(frame.faces.empty())
            throw std::runtime_error("pump relay selected no exterior faces");
        double totalArea = 0.0;
        for(auto const& face : frame.faces)
        {
            frame.origin = frame.origin + face.centroid * face.area;
            frame.normal = frame.normal + face.normal * face.area;
            totalArea += face.area;
        }
        frame.origin = frame.origin * (1.0 / totalArea);
        frame.normal = hostNormalize(frame.normal);
        frame.u = hostNormalize(
            (frame.faces.front().vertices[1] - frame.faces.front().vertices[0])
            - frame.normal
                  * hase::core::dot(frame.faces.front().vertices[1] - frame.faces.front().vertices[0], frame.normal));
        if(frame.u.euclidLength() == 0.0)
            frame.u = perpendicular(frame.normal);
        frame.v = hase::core::cross(frame.normal, frame.u);
        double scale = 0.0;
        for(auto const& face : frame.faces)
            for(auto const& vertex : face.vertices)
                scale = std::max(scale, (vertex - frame.origin).euclidLength());
        for(auto const& face : frame.faces)
        {
            if(std::abs(hase::core::dot(face.centroid - frame.origin, frame.normal)) > 1.0e-8 * std::max(1.0, scale))
                throw std::runtime_error("pump relay surfaces must be coplanar");
        }
        return frame;
    }

    struct PumpRelayGeometry
    {
        std::vector<PumpRelayDeviceDescriptor> descriptors;
        std::vector<unsigned> exitMask;
        std::vector<unsigned> entryFaceIds;
        unsigned faceCount = 0u;
    };

    [[nodiscard]] inline PumpRelayGeometry preparePumpRelayGeometry(
        hase::core::HostMesh const& mesh,
        std::vector<hase::core::PumpRelayParameters> const& relays)
    {
        PumpRelayGeometry result;
        result.faceCount = mesh.numberOfCells * mesh.numberOfFacesPerCell;
        result.exitMask.assign(relays.size() * result.faceCount, 0u);
        result.descriptors.reserve(relays.size());
        for(std::size_t relayIndex = 0u; relayIndex < relays.size(); ++relayIndex)
        {
            auto const& relay = relays[relayIndex];
            auto const exitFrame = makeRelayFrame(mesh, relay.exitSurfaces);
            auto const entryFrame = makeRelayFrame(mesh, relay.entrySurfaces);
            PumpRelayDeviceDescriptor descriptor;
            descriptor.exitOrigin = exitFrame.origin;
            descriptor.exitU = exitFrame.u;
            descriptor.exitV = exitFrame.v;
            descriptor.exitNormal = exitFrame.normal;
            descriptor.entryOrigin = entryFrame.origin;
            descriptor.entryU = entryFrame.u;
            descriptor.entryV = entryFrame.v;
            descriptor.entryNormal = entryFrame.normal;
            descriptor.cosine = std::cos(relay.rotation);
            descriptor.sine = std::sin(relay.rotation);
            descriptor.offsetU = relay.offset[0u];
            descriptor.offsetV = relay.offset[1u];
            descriptor.tiltU = relay.tilt[0u];
            descriptor.tiltV = relay.tilt[1u];
            descriptor.magnification = relay.magnification;
            descriptor.transmission = relay.transmission;
            descriptor.flipU = relay.flipU ? -1 : 1;
            descriptor.flipV = relay.flipV ? -1 : 1;
            descriptor.entryFaceBegin = static_cast<unsigned>(result.entryFaceIds.size());
            for(auto const& face : entryFrame.faces)
                result.entryFaceIds.push_back(face.cell * mesh.numberOfFacesPerCell + face.localFace);
            descriptor.entryFaceEnd = static_cast<unsigned>(result.entryFaceIds.size());
            for(auto const& face : exitFrame.faces)
            {
                unsigned const faceId = face.cell * mesh.numberOfFacesPerCell + face.localFace;
                result.exitMask[relayIndex * result.faceCount + faceId] = 1u;
            }
            result.descriptors.push_back(descriptor);
        }
        return result;
    }

    template<typename T_Value, typename T_Getter>
    [[nodiscard]] inline std::vector<T_Value> pumpRayValues(
        std::vector<GeneralPumpRayState> const& rays,
        T_Getter getter)
    {
        std::vector<T_Value> result;
        result.reserve(rays.size());
        for(auto const& ray : rays)
            result.push_back(static_cast<T_Value>(getter(ray)));
        if(result.empty())
            result.emplace_back();
        return result;
    }

    template<typename T_Value>
    [[nodiscard]] inline std::vector<T_Value> pumpDeviceStorage(std::vector<T_Value> values)
    {
        if(values.empty())
            values.emplace_back();
        return values;
    }

    template<typename T_Device>
    class GeneralPumpDeviceSource
    {
        using T_DoubleBuffer
            = ALPAKA_TYPEOF(alpaka::onHost::alloc<double>(std::declval<T_Device&>(), std::size_t{1u}));
        using T_UnsignedBuffer
            = ALPAKA_TYPEOF(alpaka::onHost::alloc<unsigned>(std::declval<T_Device&>(), std::size_t{1u}));
        using T_IntBuffer = ALPAKA_TYPEOF(alpaka::onHost::alloc<int>(std::declval<T_Device&>(), std::size_t{1u}));
        using T_DescriptorBuffer = ALPAKA_TYPEOF(
            alpaka::onHost::alloc<PumpRelayDeviceDescriptor>(std::declval<T_Device&>(), std::size_t{1u}));

    public:
        GeneralPumpDeviceSource(
            auto const& queue,
            std::vector<GeneralPumpRayState> const& rays,
            PumpRelayGeometry geometry,
            unsigned const pumpSteps)
            : m_originX(
                  hase::alpakaUtils::toDevice(
                      queue,
                      pumpRayValues<double>(rays, [](auto const& ray) { return ray.position.x; })))
            , m_originY(
                  hase::alpakaUtils::toDevice(
                      queue,
                      pumpRayValues<double>(rays, [](auto const& ray) { return ray.position.y; })))
            , m_originZ(
                  hase::alpakaUtils::toDevice(
                      queue,
                      pumpRayValues<double>(rays, [](auto const& ray) { return ray.position.z; })))
            , m_directionX(
                  hase::alpakaUtils::toDevice(
                      queue,
                      pumpRayValues<double>(rays, [](auto const& ray) { return ray.direction.x; })))
            , m_directionY(
                  hase::alpakaUtils::toDevice(
                      queue,
                      pumpRayValues<double>(rays, [](auto const& ray) { return ray.direction.y; })))
            , m_directionZ(
                  hase::alpakaUtils::toDevice(
                      queue,
                      pumpRayValues<double>(rays, [](auto const& ray) { return ray.direction.z; })))
            , m_power(
                  hase::alpakaUtils::toDevice(
                      queue,
                      pumpRayValues<double>(rays, [](auto const& ray) { return ray.power; })))
            , m_wavelength(
                  hase::alpakaUtils::toDevice(
                      queue,
                      pumpRayValues<double>(rays, [](auto const& ray) { return ray.wavelength; })))
            , m_sigmaAbsorption(
                  hase::alpakaUtils::toDevice(
                      queue,
                      pumpRayValues<double>(rays, [](auto const& ray) { return ray.sigmaAbsorption; })))
            , m_sigmaEmission(
                  hase::alpakaUtils::toDevice(
                      queue,
                      pumpRayValues<double>(rays, [](auto const& ray) { return ray.sigmaEmission; })))
            , m_cell(
                  hase::alpakaUtils::toDevice(
                      queue,
                      pumpRayValues<unsigned>(rays, [](auto const& ray) { return ray.cell; })))
            , m_forbiddenFace(
                  hase::alpakaUtils::toDevice(
                      queue,
                      pumpRayValues<int>(rays, [](auto const& ray) { return ray.forbiddenFace; })))
            , m_descriptors(hase::alpakaUtils::toDevice(queue, pumpDeviceStorage(geometry.descriptors)))
            , m_exitMask(hase::alpakaUtils::toDevice(queue, pumpDeviceStorage(geometry.exitMask)))
            , m_entryFaceIds(hase::alpakaUtils::toDevice(queue, pumpDeviceStorage(geometry.entryFaceIds)))
            , m_cacheState(
                  hase::alpakaUtils::toDevice(
                      queue,
                      std::vector<unsigned>(std::max<std::size_t>(1u, geometry.descriptors.size() * rays.size()), 0u)))
            , m_cacheTargetFace(
                  hase::alpakaUtils::toDevice(
                      queue,
                      std::vector<unsigned>(std::max<std::size_t>(1u, geometry.descriptors.size() * rays.size()), 0u)))
            , m_cacheBarycentric0(
                  hase::alpakaUtils::toDevice(
                      queue,
                      std::vector<double>(std::max<std::size_t>(1u, geometry.descriptors.size() * rays.size()), 0.0)))
            , m_cacheBarycentric1(
                  hase::alpakaUtils::toDevice(
                      queue,
                      std::vector<double>(std::max<std::size_t>(1u, geometry.descriptors.size() * rays.size()), 0.0)))
            , m_cacheBarycentric2(
                  hase::alpakaUtils::toDevice(
                      queue,
                      std::vector<double>(std::max<std::size_t>(1u, geometry.descriptors.size() * rays.size()), 0.0)))
            , m_faceCount{geometry.faceCount}
            , m_relayCount{static_cast<unsigned>(geometry.descriptors.size())}
            , m_rayCount{static_cast<unsigned>(rays.size())}
            , m_pumpSteps{pumpSteps}
        {
        }

        [[nodiscard]] unsigned rayCount() const
        {
            return m_rayCount;
        }

        [[nodiscard]] bool active(unsigned const simulationStep) const
        {
            return simulationStep < m_pumpSteps;
        }

        [[nodiscard]] auto const& cacheState() const
        {
            return m_cacheState;
        }

        void enqueue(
            auto& devBundle,
            auto const& queue,
            hase::core::DeviceMeshView const mesh,
            auto& betaVolume,
            auto& vertexPumpIntegral,
            pumpBoundaryPolicy::Relay)
        {
            enqueueWithFactory(
                devBundle,
                queue,
                mesh,
                betaVolume,
                vertexPumpIntegral,
                DevicePumpRelayBoundaryFactory{
                    m_descriptors.getView(),
                    m_exitMask.getView(),
                    m_entryFaceIds.getView(),
                    m_cacheState.getView(),
                    m_cacheTargetFace.getView(),
                    m_cacheBarycentric0.getView(),
                    m_cacheBarycentric1.getView(),
                    m_cacheBarycentric2.getView(),
                    m_faceCount,
                    m_relayCount,
                    m_rayCount});
        }

        void enqueue(
            auto& devBundle,
            auto const& queue,
            hase::core::DeviceMeshView const mesh,
            auto& betaVolume,
            auto& vertexPumpIntegral,
            pumpBoundaryPolicy::SrmBarycentric)
        {
            enqueueWithFactory(
                devBundle,
                queue,
                mesh,
                betaVolume,
                vertexPumpIntegral,
                PumpBarycentricSrmBoundaryFactory{
                    m_descriptors.getView(),
                    m_exitMask.getView(),
                    m_cacheState.getView(),
                    m_cacheTargetFace.getView(),
                    m_cacheBarycentric0.getView(),
                    m_cacheBarycentric1.getView(),
                    m_cacheBarycentric2.getView(),
                    m_faceCount,
                    m_relayCount,
                    m_rayCount});
        }

    private:
        void enqueueWithFactory(
            auto& devBundle,
            auto const& queue,
            hase::core::DeviceMeshView const mesh,
            auto& betaVolume,
            auto& vertexPumpIntegral,
            auto boundaryPolicyFactory)
        {
            if(m_rayCount == 0u)
                return;
            auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                devBundle.device,
                devBundle.executor,
                alpaka::Vec{m_rayCount});
            queue.enqueue(
                frameSpec,
                alpaka::KernelBundle{
                    TraceGeneralPump{},
                    mesh,
                    betaVolume,
                    m_originX,
                    m_originY,
                    m_originZ,
                    m_directionX,
                    m_directionY,
                    m_directionZ,
                    m_power,
                    m_wavelength,
                    m_sigmaAbsorption,
                    m_sigmaEmission,
                    m_cell,
                    m_forbiddenFace,
                    boundaryPolicyFactory,
                    vertexPumpIntegral,
                    m_rayCount});
        }

        T_DoubleBuffer m_originX, m_originY, m_originZ;
        T_DoubleBuffer m_directionX, m_directionY, m_directionZ;
        T_DoubleBuffer m_power, m_wavelength, m_sigmaAbsorption, m_sigmaEmission;
        T_UnsignedBuffer m_cell;
        T_IntBuffer m_forbiddenFace;
        T_DescriptorBuffer m_descriptors;
        T_UnsignedBuffer m_exitMask, m_entryFaceIds, m_cacheState, m_cacheTargetFace;
        T_DoubleBuffer m_cacheBarycentric0, m_cacheBarycentric1, m_cacheBarycentric2;
        unsigned m_faceCount = 0u;
        unsigned m_relayCount = 0u;
        unsigned m_rayCount = 0u;
        unsigned m_pumpSteps = 0u;
    };

    template<typename T_Device>
    [[nodiscard]] inline std::vector<GeneralPumpDeviceSource<T_Device>> prepareGeneralPumpDeviceSources(
        auto const& queue,
        hase::core::HostMesh const& mesh,
        hase::core::PumpParameters const& pump,
        unsigned const firstRay = 0u,
        unsigned const localRayCount = std::numeric_limits<unsigned>::max())
    {
        std::vector<GeneralPumpDeviceSource<T_Device>> result;
        result.reserve(pump.sources.size());
        for(auto const& source : pump.sources)
        {
            if(source.pumpSteps == 0u)
            {
                result.emplace_back(queue, std::vector<GeneralPumpRayState>{}, PumpRelayGeometry{}, 0u);
                continue;
            }
            auto rays = samplePumpSource(mesh, source, source.rayCount, source.rngSeed, firstRay, localRayCount);
            result.emplace_back(queue, rays, preparePumpRelayGeometry(mesh, source.relays), source.pumpSteps);
        }
        return result;
    }

    struct ProjectVertexPumpRateToCells
    {
        ALPAKA_FN_ACC void operator()(
            auto const& acc,
            hase::core::DeviceMeshView const mesh,
            auto vertexIntegral,
            auto lumpedVertexVolume,
            auto cellRate) const
        {
            for(auto [cell] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfCells}))
            {
                if(mesh.getCellType(cell) == mesh.claddingNumber)
                {
                    cellRate[cell] = 0.0;
                    continue;
                }

                double rateSum = 0.0;
                for(unsigned localVertex = 0u; localVertex < mesh.numberOfCellVertices; ++localVertex)
                {
                    unsigned const point = mesh.cellPointIndices[cell * mesh.numberOfCellVertices + localVertex];
                    double const volume = lumpedVertexVolume[point];
                    rateSum += volume > 0.0 ? vertexIntegral[point] / volume : 0.0;
                }
                // Volume averaging the four lumped vertex rates preserves the deposited
                // integral: sum_cell(V_cell * cellRate) == sum_vertex(vertexIntegral).
                cellRate[cell] = rateSum / static_cast<double>(mesh.numberOfCellVertices);
            }
        }
    };

    template<
        typename T_Device,
        typename T_Executor,
        typename T_BetaBuffer,
        typename T_VertexBuffer,
        typename T_BoundaryPolicy = pumpBoundaryPolicy::Relay>
    void enqueueGeneralPumpIntegrals(
        hase::alpakaUtils::DevBundle<T_Device, T_Executor>& devBundle,
        auto const& queue,
        hase::core::DeviceMeshView const mesh,
        std::vector<GeneralPumpDeviceSource<T_Device>>& sources,
        T_BetaBuffer& betaVolume,
        T_VertexBuffer& vertexPumpIntegral,
        unsigned const simulationStep,
        T_BoundaryPolicy boundaryPolicy = {})
    {
        HASE_HOST_ROUTINE_SCOPE("pump.enqueue_integrals");
        alpaka::onHost::fill(
            queue,
            vertexPumpIntegral,
            0.0,
            alpaka::Vec{static_cast<std::size_t>(mesh.numberOfMeshPoints)});
        for(auto& source : sources)
        {
            if(source.active(simulationStep))
                source.enqueue(devBundle, queue, mesh, betaVolume, vertexPumpIntegral, boundaryPolicy);
        }
    }

    template<
        typename T_Device,
        typename T_Executor,
        typename T_VertexBuffer,
        typename T_LumpedVolumeBuffer,
        typename T_RateBuffer>
    void enqueueProjectVertexPumpRateToCells(
        hase::alpakaUtils::DevBundle<T_Device, T_Executor>& devBundle,
        auto const& queue,
        hase::core::DeviceMeshView const mesh,
        T_VertexBuffer& vertexPumpIntegral,
        T_LumpedVolumeBuffer& lumpedVertexVolume,
        T_RateBuffer& cellRate)
    {
        HASE_HOST_ROUTINE_SCOPE("pump.enqueue_project_vertices");
        auto cellFrameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
            devBundle.device,
            devBundle.executor,
            alpaka::Vec{mesh.numberOfCells});
        queue.enqueue(
            cellFrameSpec,
            alpaka::KernelBundle{
                ProjectVertexPumpRateToCells{},
                mesh,
                vertexPumpIntegral,
                lumpedVertexVolume,
                cellRate});
    }

    template<
        typename T_Device,
        typename T_Executor,
        typename T_BetaBuffer,
        typename T_VertexBuffer,
        typename T_LumpedVolumeBuffer,
        typename T_RateBuffer,
        typename T_BoundaryPolicy = pumpBoundaryPolicy::Relay>
    void enqueueGeneralPump(
        hase::alpakaUtils::DevBundle<T_Device, T_Executor>& devBundle,
        auto const& queue,
        hase::core::DeviceMeshView const mesh,
        std::vector<GeneralPumpDeviceSource<T_Device>>& sources,
        T_BetaBuffer& betaVolume,
        T_VertexBuffer& vertexPumpIntegral,
        T_LumpedVolumeBuffer& lumpedVertexVolume,
        T_RateBuffer& cellRate,
        unsigned const simulationStep,
        T_BoundaryPolicy boundaryPolicy = {})
    {
        enqueueGeneralPumpIntegrals(
            devBundle,
            queue,
            mesh,
            sources,
            betaVolume,
            vertexPumpIntegral,
            simulationStep,
            boundaryPolicy);
        enqueueProjectVertexPumpRateToCells(devBundle, queue, mesh, vertexPumpIntegral, lumpedVertexVolume, cellRate);
    }
} // namespace hase::kernels
