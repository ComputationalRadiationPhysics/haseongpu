/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <kernels/forward/rayWalk.hpp>

#include <cstdint>
#include <limits>

namespace hase::kernels::forward
{
    inline constexpr unsigned maxImmediateFaceTransitions = 32u;
    inline constexpr double ownershipProbeEpsilon = 1.0e-12;

    enum class Tet4TransitionStatus : std::uint8_t
    {
        enteredCell,
        reachedBoundary,
        failed
    };

    struct Tet4FaceTransition
    {
        unsigned cell = 0u;
        int forbiddenFace = -1;
        int boundaryFace = -1;
        Tet4TransitionStatus status = Tet4TransitionStatus::enteredCell;
    };

    [[nodiscard]] inline ALPAKA_FN_ACC bool isNearTet4Face(
        hase::core::DeviceMeshView const& mesh,
        unsigned const tet,
        hase::core::Point const point)
    {
        for(unsigned localFace = 0u; localFace < mesh.numberOfFacesPerCell; ++localFace)
        {
            if(alpaka::math::abs(mesh.getFaceBarycentricCoordinate(tet, localFace, point))
               <= barycentricTraversalTolerance)
            {
                return true;
            }
        }
        return false;
    }

    [[nodiscard]] inline ALPAKA_FN_ACC int immediateExitFace(
        hase::core::DeviceMeshView const& mesh,
        unsigned const tet,
        hase::core::Point const point,
        hase::core::Point const direction,
        int const forbiddenFace)
    {
        for(unsigned localFace = 0u; localFace < mesh.numberOfFacesPerCell; ++localFace)
        {
            if(static_cast<int>(localFace) == forbiddenFace)
            {
                continue;
            }
            double const coordinate = mesh.getFaceBarycentricCoordinate(tet, localFace, point);
            double const directionalChange = mesh.getFaceBarycentricDirection(tet, localFace, direction);
            if(alpaka::math::abs(coordinate) <= barycentricTraversalTolerance && directionalChange < 0.0)
            {
                return static_cast<int>(localFace);
            }
        }
        return -1;
    }

    [[nodiscard]] inline ALPAKA_FN_ACC Tet4FaceTransition recoverFaceTransition(
        hase::core::DeviceMeshView const& mesh,
        unsigned cell,
        int exitFace,
        hase::core::Point const point,
        hase::core::Point const direction)
    {
        Tet4FaceTransition result{cell};
        constexpr unsigned invalidCell = std::numeric_limits<unsigned>::max();
        unsigned recentCell0 = cell;
        unsigned recentCell1 = invalidCell;
        unsigned recentCell2 = invalidCell;
        for(unsigned transition = 0u; transition < maxImmediateFaceTransitions; ++transition)
        {
            int const neighbor = mesh.getCellNeighbor(cell, static_cast<unsigned>(exitFace));
            if(neighbor < 0)
            {
                result.cell = cell;
                result.boundaryFace = exitFace;
                result.status = Tet4TransitionStatus::reachedBoundary;
                return result;
            }
            int const forbiddenFace = mesh.getCellNeighborLocalFace(cell, static_cast<unsigned>(exitFace));
            if(static_cast<unsigned>(neighbor) >= mesh.numberOfCells || forbiddenFace < 0
               || static_cast<unsigned>(forbiddenFace) >= mesh.numberOfFacesPerCell)
            {
                result.cell = cell;
                result.status = Tet4TransitionStatus::failed;
                return result;
            }

            cell = static_cast<unsigned>(neighbor);
            result.cell = cell;
            result.forbiddenFace = forbiddenFace;

            int fallbackExitFace = -1;
            int preferredExitFace = -1;
            recentCell2 = recentCell1;
            recentCell1 = recentCell0;
            recentCell0 = cell;
            for(unsigned localFace = 0u; localFace < mesh.numberOfFacesPerCell; ++localFace)
            {
                // Test point + epsilon * direction in barycentric space without moving the ray.
                double coordinate = mesh.getFaceBarycentricCoordinate(cell, localFace, point);
                if(alpaka::math::abs(coordinate) <= barycentricTraversalTolerance)
                {
                    coordinate = 0.0;
                }
                double const directionalChange = mesh.getFaceBarycentricDirection(cell, localFace, direction);
                double const probeCoordinate = coordinate + ownershipProbeEpsilon * directionalChange;
                if(probeCoordinate >= 0.0)
                {
                    continue;
                }

                if(fallbackExitFace < 0)
                {
                    fallbackExitFace = static_cast<int>(localFace);
                }
                int const candidateNeighbor = mesh.getCellNeighbor(cell, localFace);
                if(candidateNeighbor < 0)
                {
                    continue;
                }
                unsigned const candidateCell = static_cast<unsigned>(candidateNeighbor);
                if(candidateCell != recentCell0 && candidateCell != recentCell1 && candidateCell != recentCell2)
                {
                    preferredExitFace = static_cast<int>(localFace);
                    break;
                }
            }

            exitFace = preferredExitFace >= 0 ? preferredExitFace : fallbackExitFace;
            if(exitFace < 0)
            {
                return result;
            }
        }
        result.status = Tet4TransitionStatus::failed;
        return result;
    }

    [[nodiscard]] inline ALPAKA_FN_ACC Tet4FaceTransition transitionAcrossIntersection(
        hase::core::DeviceMeshView const& mesh,
        unsigned const cell,
        Tet4FaceIntersection const intersection,
        hase::core::Point const point,
        hase::core::Point const direction)
    {
        if(intersection.localFace < 0)
        {
            return Tet4FaceTransition{cell, -1, -1, Tet4TransitionStatus::failed};
        }
        if(hasMultipleTiedFaces(intersection.tiedFaceMask))
        {
            return recoverFaceTransition(mesh, cell, intersection.localFace, point, direction);
        }

        int const neighbor = mesh.getCellNeighbor(cell, static_cast<unsigned>(intersection.localFace));
        if(neighbor < 0)
        {
            return Tet4FaceTransition{cell, -1, intersection.localFace, Tet4TransitionStatus::reachedBoundary};
        }
        int const forbiddenFace = mesh.getCellNeighborLocalFace(cell, static_cast<unsigned>(intersection.localFace));
        if(static_cast<unsigned>(neighbor) >= mesh.numberOfCells || forbiddenFace < 0
           || static_cast<unsigned>(forbiddenFace) >= mesh.numberOfFacesPerCell)
        {
            return Tet4FaceTransition{cell, -1, -1, Tet4TransitionStatus::failed};
        }
        return Tet4FaceTransition{
            static_cast<unsigned>(neighbor),
            forbiddenFace,
            -1,
            Tet4TransitionStatus::enteredCell};
    }
} // namespace hase::kernels::forward
