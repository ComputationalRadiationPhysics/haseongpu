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
#include <core/mesh.hpp>
#include <core/simulationRunControl.hpp>
#include <kernels/forward/rayTransition.hpp>
#include <kernels/forward/rayWalk.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
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

    struct PumpRayBatch
    {
        std::vector<double> originX, originY, originZ;
        std::vector<double> directionX, directionY, directionZ;
        std::vector<double> power, wavelength, sigmaAbsorption, sigmaEmission;
        std::vector<unsigned> cell;
        std::vector<int> forbiddenFace, exitFace;

        [[nodiscard]] std::size_t size() const
        {
            return power.size();
        }
    };

    [[nodiscard]] inline bool containsDomain(std::vector<int> const& domains, int const domain)
    {
        return std::ranges::find(domains, domain) != domains.end();
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

    [[nodiscard]] inline PumpRayBatch samplePumpSource(
        hase::core::HostMesh const& mesh,
        hase::core::PumpSourceParameters const& source,
        unsigned const rayCount,
        std::uint32_t const seed)
    {
        auto const faces = pumpBoundaryFaces(mesh, source.surfaces);
        if(faces.empty())
            throw std::runtime_error("pump source selected no exterior boundary faces");
        std::vector<double> areas;
        areas.reserve(faces.size());
        for(auto const& face : faces)
            areas.push_back(face.area);
        std::discrete_distribution<std::size_t> faceDistribution(areas.begin(), areas.end());
        std::discrete_distribution<std::size_t> spectrumDistribution(
            source.spectralWeights.begin(),
            source.spectralWeights.end());
        std::discrete_distribution<std::size_t> angularDistribution(
            source.angularWeights.begin(),
            source.angularWeights.end());
        std::mt19937_64 rng(seed);
        std::uniform_real_distribution<double> uniform(0.0, 1.0);

        PumpRayBatch batch;
        auto reserve = [rayCount](auto& values) { values.reserve(rayCount); };
        reserve(batch.originX);
        reserve(batch.originY);
        reserve(batch.originZ);
        reserve(batch.directionX);
        reserve(batch.directionY);
        reserve(batch.directionZ);
        reserve(batch.power);
        reserve(batch.wavelength);
        reserve(batch.sigmaAbsorption);
        reserve(batch.sigmaEmission);
        reserve(batch.cell);
        reserve(batch.forbiddenFace);
        reserve(batch.exitFace);
        for(unsigned ray = 0u; ray < rayCount; ++ray)
        {
            PumpBoundaryFace const* face = nullptr;
            hase::core::Point origin;
            bool accepted = false;
            for(unsigned attempt = 0u; attempt < 100000u; ++attempt)
            {
                face = &faces[faceDistribution(rng)];
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

            batch.originX.push_back(origin.x);
            batch.originY.push_back(origin.y);
            batch.originZ.push_back(origin.z);
            batch.directionX.push_back(direction.x);
            batch.directionY.push_back(direction.y);
            batch.directionZ.push_back(direction.z);
            batch.power.push_back(source.totalPower / static_cast<double>(rayCount));
            batch.wavelength.push_back(source.wavelengths[spectrum]);
            batch.sigmaAbsorption.push_back(source.sigmaAbsorption[spectrum]);
            batch.sigmaEmission.push_back(source.sigmaEmission[spectrum]);
            batch.cell.push_back(face->cell);
            batch.forbiddenFace.push_back(static_cast<int>(face->localFace));
            batch.exitFace.push_back(-1);
        }
        return batch;
    }

    struct TraceGeneralPump
    {
        double planckConstant = 6.62607015e-34;
        double speedOfLight = 299792458.0;

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
            typename T_ExitFaceView,
            typename T_CellPumpIntegralView,
            typename T_SamplePumpIntegralView>
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
            T_ExitFaceView exitFace,
            T_CellPumpIntegralView cellPumpIntegral,
            T_SamplePumpIntegralView samplePumpIntegral,
            unsigned const rayCount) const
        {
            for(auto [ray] :
                alpaka::onAcc::makeIdxMap(acc, alpaka::onAcc::worker::threadsInGrid, alpaka::IdxRange{rayCount}))
            {
                hase::core::Point origin{originX[ray], originY[ray], originZ[ray]};
                hase::core::Point const direction{directionX[ray], directionY[ray], directionZ[ray]};
                unsigned tet = cell[ray];
                int forbidden = forbiddenFace[ray];
                double rayPower = power[ray];
                exitFace[ray] = -1;
                constexpr unsigned maxSteps = 10000u;
                for(unsigned step = 0u; step < maxSteps && rayPower != 0.0; ++step)
                {
                    auto const intersection
                        = hase::kernels::forward::nextFaceIntersection(mesh, tet, origin, direction, forbidden);
                    if(intersection.localFace < 0)
                    {
                        rayPower = 0.0;
                        break;
                    }
                    bool const gainCell = mesh.getCellType(tet) != mesh.claddingNumber;
                    double const gain = gainCell ? static_cast<double>(mesh.nTot)
                                                       * (betaVolume[tet] * (sigmaAbsorption[ray] + sigmaEmission[ray])
                                                          - sigmaAbsorption[ray])
                                                 : -mesh.claddingAbsorption;
                    double const exponent = gain * intersection.length;
                    if(!alpaka::math::isfinite(exponent) || exponent > 700.0)
                    {
                        rayPower = 0.0;
                        break;
                    }
                    double const nextPower = rayPower * alpaka::math::exp(exponent);
                    if(gainCell && mesh.nTot > 0.0f)
                    {
                        double const integral = (rayPower - nextPower) * wavelength[ray]
                                                / (planckConstant * speedOfLight * static_cast<double>(mesh.nTot));
                        alpaka::onAcc::atomicAdd(acc, &cellPumpIntegral[tet], integral);
                        if(mesh.samplePointsAreMeshPoints)
                        {
                            auto const midpoint = origin + direction * (0.5 * intersection.length);
                            auto const barycentric
                                = hase::kernels::forward::barycentricCoordinates(mesh, tet, midpoint);
                            for(unsigned vertex = 0u; vertex < mesh.numberOfCellVertices; ++vertex)
                                alpaka::onAcc::atomicAdd(
                                    acc,
                                    &samplePumpIntegral
                                        [mesh.cellPointIndices[tet * mesh.numberOfCellVertices + vertex]],
                                    integral * barycentric[vertex]);
                        }
                    }
                    rayPower = nextPower;
                    origin = hase::kernels::forward::advance(origin, direction, intersection.length);
                    auto const transition = hase::kernels::forward::transitionAcrossIntersection(
                        mesh,
                        tet,
                        intersection,
                        origin,
                        direction);
                    if(transition.status == hase::kernels::forward::Tet4TransitionStatus::failed)
                    {
                        rayPower = 0.0;
                        break;
                    }
                    if(transition.status == hase::kernels::forward::Tet4TransitionStatus::reachedBoundary)
                    {
                        exitFace[ray] = transition.boundaryFace;
                        break;
                    }
                    tet = transition.cell;
                    forbidden = transition.forbiddenFace;
                }
                originX[ray] = origin.x;
                originY[ray] = origin.y;
                originZ[ray] = origin.z;
                power[ray] = rayPower;
                cell[ray] = tet;
                forbiddenFace[ray] = forbidden;
            }
        }
    };

    template<
        typename T_Device,
        typename T_Executor,
        typename T_BetaBuffer,
        typename T_CellBuffer,
        typename T_SampleBuffer>
    PumpRayBatch tracePumpBatch(
        hase::alpakaUtils::DevBundle<T_Device, T_Executor>& devBundle,
        auto const& queue,
        hase::core::DeviceMeshView const mesh,
        T_BetaBuffer& betaVolume,
        T_CellBuffer& cellPumpIntegral,
        T_SampleBuffer& samplePumpIntegral,
        PumpRayBatch batch)
    {
        unsigned const count = static_cast<unsigned>(batch.size());
        if(count == 0u)
            return batch;
        auto originX = hase::alpakaUtils::toDevice(queue, batch.originX);
        auto originY = hase::alpakaUtils::toDevice(queue, batch.originY);
        auto originZ = hase::alpakaUtils::toDevice(queue, batch.originZ);
        auto directionX = hase::alpakaUtils::toDevice(queue, batch.directionX);
        auto directionY = hase::alpakaUtils::toDevice(queue, batch.directionY);
        auto directionZ = hase::alpakaUtils::toDevice(queue, batch.directionZ);
        auto power = hase::alpakaUtils::toDevice(queue, batch.power);
        auto wavelength = hase::alpakaUtils::toDevice(queue, batch.wavelength);
        auto sigmaA = hase::alpakaUtils::toDevice(queue, batch.sigmaAbsorption);
        auto sigmaE = hase::alpakaUtils::toDevice(queue, batch.sigmaEmission);
        auto cell = hase::alpakaUtils::toDevice(queue, batch.cell);
        auto forbiddenFace = hase::alpakaUtils::toDevice(queue, batch.forbiddenFace);
        auto exitFace = hase::alpakaUtils::toDevice(queue, batch.exitFace);
        auto frameSpec
            = hase::alpakaUtils::getFrameSpec<uint32_t>(devBundle.device, devBundle.executor, alpaka::Vec{count});
        queue.enqueue(
            frameSpec,
            alpaka::KernelBundle{
                TraceGeneralPump{},
                mesh,
                betaVolume,
                originX,
                originY,
                originZ,
                directionX,
                directionY,
                directionZ,
                power,
                wavelength,
                sigmaA,
                sigmaE,
                cell,
                forbiddenFace,
                exitFace,
                cellPumpIntegral,
                samplePumpIntegral,
                count});
        alpaka::onHost::wait(queue);
        auto copyBack = [&](auto const& deviceBuffer, auto& values)
        {
            auto host = alpaka::onHost::allocHostLike(deviceBuffer);
            alpaka::onHost::memcpy(queue, host, deviceBuffer);
            alpaka::onHost::wait(queue);
            std::copy_n(alpaka::onHost::data(host), values.size(), values.begin());
        };
        copyBack(originX, batch.originX);
        copyBack(originY, batch.originY);
        copyBack(originZ, batch.originZ);
        copyBack(directionX, batch.directionX);
        copyBack(directionY, batch.directionY);
        copyBack(directionZ, batch.directionZ);
        copyBack(power, batch.power);
        copyBack(cell, batch.cell);
        copyBack(forbiddenFace, batch.forbiddenFace);
        copyBack(exitFace, batch.exitFace);
        return batch;
    }

    struct RelayFrame
    {
        hase::core::Point origin, u, v, normal;
        std::vector<PumpBoundaryFace> faces;
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

    [[nodiscard]] inline bool pointInTriangle(
        hase::core::Point const point,
        PumpBoundaryFace const& face,
        hase::core::Point const u,
        hase::core::Point const v)
    {
        auto project = [&](hase::core::Point const p)
        { return std::array<double, 2u>{hase::core::dot(p, u), hase::core::dot(p, v)}; };
        auto const p = project(point);
        auto const a = project(face.vertices[0]);
        auto const b = project(face.vertices[1]);
        auto const c = project(face.vertices[2]);
        double const denominator = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1]);
        if(std::abs(denominator) <= 1.0e-30)
            return false;
        double const wa = ((b[1] - c[1]) * (p[0] - c[0]) + (c[0] - b[0]) * (p[1] - c[1])) / denominator;
        double const wb = ((c[1] - a[1]) * (p[0] - c[0]) + (a[0] - c[0]) * (p[1] - c[1])) / denominator;
        double const wc = 1.0 - wa - wb;
        return wa >= -1.0e-10 && wb >= -1.0e-10 && wc >= -1.0e-10;
    }

    [[nodiscard]] inline PumpRayBatch applyPumpRelay(
        hase::core::HostMesh const& mesh,
        PumpRayBatch const& exits,
        hase::core::PumpRelayParameters const& relay)
    {
        auto const exitFrame = makeRelayFrame(mesh, relay.exitSurfaces);
        auto const entryFrame = makeRelayFrame(mesh, relay.entrySurfaces);
        PumpRayBatch result;
        double const cosine = std::cos(relay.rotation);
        double const sine = std::sin(relay.rotation);
        for(std::size_t ray = 0u; ray < exits.size(); ++ray)
        {
            if(exits.exitFace[ray] < 0 || exits.power[ray] == 0.0)
                continue;
            unsigned const faceIndex
                = exits.cell[ray] * mesh.numberOfFacesPerCell + static_cast<unsigned>(exits.exitFace[ray]);
            int const domain = mesh.cellFaceBoundaries[faceIndex];
            if(!containsDomain(relay.exitSurfaces, domain))
                continue;
            hase::core::Point const position{exits.originX[ray], exits.originY[ray], exits.originZ[ray]};
            hase::core::Point const relative = position - exitFrame.origin;
            double u = hase::core::dot(relative, exitFrame.u) * (relay.flipU ? -1.0 : 1.0);
            double v = hase::core::dot(relative, exitFrame.v) * (relay.flipV ? -1.0 : 1.0);
            u *= relay.magnification;
            v *= relay.magnification;
            double const mappedU = cosine * u - sine * v + relay.offset[0];
            double const mappedV = sine * u + cosine * v + relay.offset[1];
            hase::core::Point const mappedPosition
                = entryFrame.origin + entryFrame.u * mappedU + entryFrame.v * mappedV;

            PumpBoundaryFace const* entryFace = nullptr;
            for(auto const& candidate : entryFrame.faces)
            {
                if(pointInTriangle(mappedPosition, candidate, entryFrame.u, entryFrame.v))
                {
                    entryFace = &candidate;
                    break;
                }
            }
            if(entryFace == nullptr)
                continue;

            hase::core::Point const oldDirection{exits.directionX[ray], exits.directionY[ray], exits.directionZ[ray]};
            double du = hase::core::dot(oldDirection, exitFrame.u) * (relay.flipU ? -1.0 : 1.0);
            double dv = hase::core::dot(oldDirection, exitFrame.v) * (relay.flipV ? -1.0 : 1.0);
            double const mappedDu = cosine * du - sine * dv + relay.tilt[0];
            double const mappedDv = sine * du + cosine * dv + relay.tilt[1];
            double const normalMagnitude = std::abs(hase::core::dot(oldDirection, exitFrame.normal));
            hase::core::Point const direction = hostNormalize(
                entryFrame.u * mappedDu + entryFrame.v * mappedDv - entryFrame.normal * normalMagnitude);

            result.originX.push_back(mappedPosition.x);
            result.originY.push_back(mappedPosition.y);
            result.originZ.push_back(mappedPosition.z);
            result.directionX.push_back(direction.x);
            result.directionY.push_back(direction.y);
            result.directionZ.push_back(direction.z);
            result.power.push_back(exits.power[ray] * relay.transmission);
            result.wavelength.push_back(exits.wavelength[ray]);
            result.sigmaAbsorption.push_back(exits.sigmaAbsorption[ray]);
            result.sigmaEmission.push_back(exits.sigmaEmission[ray]);
            result.cell.push_back(entryFace->cell);
            result.forbiddenFace.push_back(static_cast<int>(entryFace->localFace));
            result.exitFace.push_back(-1);
        }
        return result;
    }

    struct NormalizePumpRate
    {
        ALPAKA_FN_ACC void operator()(
            auto const& acc,
            hase::core::DeviceMeshView const mesh,
            auto cellIntegral,
            auto lumpedVolume,
            auto sampleRate) const
        {
            for(auto [sample] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfSamples}))
            {
                if(mesh.samplePointsAreMeshPoints)
                    sampleRate[sample] = lumpedVolume[sample] > 0.0 ? sampleRate[sample] / lumpedVolume[sample] : 0.0;
                else
                    sampleRate[sample] = cellIntegral[sample] / mesh.getCellVolume(sample);
            }
        }
    };

    template<
        typename T_Device,
        typename T_Executor,
        typename T_BetaBuffer,
        typename T_CellBuffer,
        typename T_LumpedBuffer,
        typename T_SampleBuffer>
    void enqueueGeneralPump(
        hase::alpakaUtils::DevBundle<T_Device, T_Executor>& devBundle,
        auto const& queue,
        hase::core::HostMesh const& hostMesh,
        hase::core::DeviceMeshView const mesh,
        hase::core::PumpParameters const& pump,
        T_BetaBuffer& betaVolume,
        T_CellBuffer& cellPumpIntegral,
        T_LumpedBuffer& lumpedVolume,
        T_SampleBuffer& sampleRate)
    {
        alpaka::onHost::fill(queue, cellPumpIntegral, 0.0, alpaka::Vec{static_cast<std::size_t>(mesh.numberOfCells)});
        alpaka::onHost::fill(queue, sampleRate, 0.0, alpaka::Vec{static_cast<std::size_t>(mesh.numberOfSamples)});
        alpaka::onHost::wait(queue);
        for(std::size_t sourceIndex = 0u; sourceIndex < pump.sources.size(); ++sourceIndex)
        {
            auto const& source = pump.sources[sourceIndex];
            PumpRayBatch rays = samplePumpSource(
                hostMesh,
                source,
                pump.rayCount,
                pump.rngSeed + static_cast<std::uint32_t>(sourceIndex));
            rays = tracePumpBatch(devBundle, queue, mesh, betaVolume, cellPumpIntegral, sampleRate, std::move(rays));
            for(auto const& relay : source.relays)
            {
                rays = applyPumpRelay(hostMesh, rays, relay);
                rays = tracePumpBatch(
                    devBundle,
                    queue,
                    mesh,
                    betaVolume,
                    cellPumpIntegral,
                    sampleRate,
                    std::move(rays));
            }
        }
        auto sampleFrameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
            devBundle.device,
            devBundle.executor,
            alpaka::Vec{mesh.numberOfSamples});
        queue.enqueue(
            sampleFrameSpec,
            alpaka::KernelBundle{NormalizePumpRate{}, mesh, cellPumpIntegral, lumpedVolume, sampleRate});
        alpaka::onHost::wait(queue);
    }
} // namespace hase::kernels
