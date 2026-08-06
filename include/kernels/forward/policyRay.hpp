/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <alpaka/alpaka.hpp>

#include <core/mesh.hpp>
#include <kernels/forward/rayTransition.hpp>

#include <array>
#include <cassert>
#include <concepts>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace hase::kernels::forward::ray
{
    using TriangleBarycentric = std::array<double, 3u>;

    [[nodiscard]] ALPAKA_FN_ACC inline TriangleBarycentric triangleBarycentricCoordinates(
        hase::core::DeviceMeshView const& mesh,
        unsigned const cell,
        unsigned const localFace,
        hase::core::Point const point)
    {
        auto const a = mesh.getPoint(static_cast<unsigned>(mesh.getCellFacePoint(cell, localFace, 0u)));
        auto const b = mesh.getPoint(static_cast<unsigned>(mesh.getCellFacePoint(cell, localFace, 1u)));
        auto const c = mesh.getPoint(static_cast<unsigned>(mesh.getCellFacePoint(cell, localFace, 2u)));
        auto const v0 = b - a;
        auto const v1 = c - a;
        auto const v2 = point - a;
        double const d00 = hase::core::dot(v0, v0);
        double const d01 = hase::core::dot(v0, v1);
        double const d11 = hase::core::dot(v1, v1);
        double const d20 = hase::core::dot(v2, v0);
        double const d21 = hase::core::dot(v2, v1);
        double const denominator = d00 * d11 - d01 * d01;
        if(alpaka::math::abs(denominator) <= std::numeric_limits<double>::epsilon())
            return {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        double const second = (d11 * d20 - d01 * d21) / denominator;
        double const third = (d00 * d21 - d01 * d20) / denominator;
        return {1.0 - second - third, second, third};
    }

    [[nodiscard]] ALPAKA_FN_ACC inline hase::core::Point positionFromTriangleBarycentric(
        hase::core::DeviceMeshView const& mesh,
        unsigned const cell,
        unsigned const localFace,
        TriangleBarycentric const coordinates)
    {
        hase::core::Point result{0.0, 0.0, 0.0};
        for(unsigned vertex = 0u; vertex < 3u; ++vertex)
            result = result
                     + mesh.getPoint(static_cast<unsigned>(mesh.getCellFacePoint(cell, localFace, vertex)))
                           * coordinates[vertex];
        return result;
    }

    namespace srmPosition
    {
        struct Centroid
        {
            static constexpr bool storesBarycentric = false;

            [[nodiscard]] ALPAKA_FN_ACC static hase::core::Point restore(
                hase::core::DeviceMeshView const& mesh,
                unsigned const cell,
                unsigned const localFace)
            {
                return faceCentroid(mesh, cell, localFace);
            }
        };

        struct Barycentric
        {
            static constexpr bool storesBarycentric = true;
        };

        inline constexpr Centroid centroid;
        inline constexpr Barycentric barycentric;
    } // namespace srmPosition

    template<typename T>
    concept SrmPositionPolicy = std::same_as<std::remove_cvref_t<T>, srmPosition::Centroid>
                                || std::same_as<std::remove_cvref_t<T>, srmPosition::Barycentric>;

    namespace behaviourDimension
    {
        struct Cell
        {
        };

        struct Boundary
        {
        };
    } // namespace behaviourDimension

    namespace trait
    {
        template<typename T>
        struct IsCellBehaviour : std::bool_constant<std::derived_from<T, behaviourDimension::Cell>>
        {
        };

        template<typename T>
        struct IsBoundaryBehaviour : std::bool_constant<std::derived_from<T, behaviourDimension::Boundary>>
        {
        };

    } // namespace trait

    template<typename T>
    inline constexpr bool isCellBehaviour_v = trait::IsCellBehaviour<std::remove_cvref_t<T>>::value;

    template<typename T>
    inline constexpr bool isBoundaryBehaviour_v = trait::IsBoundaryBehaviour<std::remove_cvref_t<T>>::value;

    namespace concepts
    {
        template<typename T>
        concept CellBehaviour = isCellBehaviour_v<T>;

        template<typename T>
        concept BoundaryBehaviour = isBoundaryBehaviour_v<T>;

        template<typename T>
        concept BehaviourTerm = CellBehaviour<T> || BoundaryBehaviour<T>;

    } // namespace concepts

    /** SRM boundary tag. Its position policy determines the ray's optional storage. */
    template<SrmPositionPolicy T_PositionPolicy = srmPosition::Centroid>
    struct BoundaryPolicySrm : behaviourDimension::Boundary
    {
        using PositionPolicy = T_PositionPolicy;
        static inline constexpr PositionPolicy positionPolicy{};
    };

    template<typename T>
    concept SrmBoundaryBehaviour
        = concepts::BoundaryBehaviour<T> && requires { typename std::remove_cvref_t<T>::PositionPolicy; };

    enum class BoundaryResult : std::uint8_t
    {
        stop,
        continueTraversal
    };

    template<typename T_Term, typename T_Acc, typename T_RayState>
    concept HasOnCell = concepts::CellBehaviour<T_Term>
                        && requires(
                            T_Term& term,
                            T_Acc const& acc,
                            hase::core::DeviceMeshView const& mesh,
                            T_RayState& rayState,
                            unsigned const cell,
                            Tet4FaceIntersection const intersection) {
                               { term(acc, mesh, rayState, cell, intersection) } -> std::same_as<bool>;
                           };

    template<typename T_Term, typename T_Acc, typename T_RayState>
    concept HasOnBoundary = concepts::BoundaryBehaviour<T_Term>
                            && requires(
                                T_Term& term,
                                T_Acc const& acc,
                                hase::core::DeviceMeshView const& mesh,
                                T_RayState& rayState,
                                unsigned const cell,
                                unsigned const localFace) {
                                   { term(acc, mesh, rayState, cell, localFace) } -> std::same_as<BoundaryResult>;
                               };

    inline constexpr BoundaryPolicySrm aseSrmPolicy{};
    inline constexpr BoundaryPolicySrm<srmPosition::Barycentric> pumpSrmPolicy{};

    struct NoSrmPositionStorage
    {
    };

    struct BarycentricSrmPositionStorage
    {
        TriangleBarycentric boundaryBarycentric{};
    };

    template<SrmPositionPolicy T_PositionPolicy>
    struct SrmPositionStorage : NoSrmPositionStorage
    {
    };

    template<>
    struct SrmPositionStorage<srmPosition::Barycentric> : BarycentricSrmPositionStorage
    {
    };

    ALPAKA_FN_ACC inline void captureSrmPosition(
        srmPosition::Centroid,
        hase::core::DeviceMeshView const&,
        unsigned,
        unsigned,
        hase::core::Point const,
        NoSrmPositionStorage&)
    {
    }

    ALPAKA_FN_ACC inline void captureSrmPosition(
        srmPosition::Barycentric,
        hase::core::DeviceMeshView const& mesh,
        unsigned const cell,
        unsigned const localFace,
        hase::core::Point const position,
        BarycentricSrmPositionStorage& storage)
    {
        storage.boundaryBarycentric = triangleBarycentricCoordinates(mesh, cell, localFace, position);
    }

    [[nodiscard]] ALPAKA_FN_ACC inline hase::core::Point restoreSrmPosition(
        srmPosition::Centroid,
        hase::core::DeviceMeshView const& mesh,
        unsigned const cell,
        unsigned const localFace)
    {
        return faceCentroid(mesh, cell, localFace);
    }

    [[nodiscard]] ALPAKA_FN_ACC inline hase::core::Point restoreSrmPosition(
        srmPosition::Barycentric,
        hase::core::DeviceMeshView const& mesh,
        unsigned const cell,
        unsigned const localFace,
        BarycentricSrmPositionStorage const& storage)
    {
        return positionFromTriangleBarycentric(mesh, cell, localFace, storage.boundaryBarycentric);
    }

    struct TraversalState
    {
        hase::core::Point position{};
        hase::core::Point direction{};
        unsigned cell = 0u;
        int forbiddenFace = -1;
    };

    template<typename T_Ray>
    concept State = requires(T_Ray rayState) {
        rayState.position;
        rayState.direction;
        rayState.cell;
        rayState.forbiddenFace;
    };

    struct BoundaryPolicyEscape : behaviourDimension::Boundary
    {
        ALPAKA_FN_ACC BoundaryResult operator()(auto const&, auto const&, auto&, unsigned, unsigned)
        {
            return BoundaryResult::stop;
        }
    };

    /**
     * Compile-time composition of ray-walk behaviour terms.
     *
     * Terms are callable objects categorized as cell or boundary behaviour.
     * Different categories can be supplied in any order; multiple terms in
     * one category are composed in the supplied order.
     */
    template<concepts::BehaviourTerm... T_Terms>
    struct RayWalkBehaviour : T_Terms...
    {
        static_assert((concepts::CellBehaviour<T_Terms> || ...), "ray walk requires cell behaviour");
        static_assert(
            (std::size_t{0u} + ... + static_cast<std::size_t>(concepts::BoundaryBehaviour<T_Terms>)) == 1u,
            "ray walk requires exactly one boundary behaviour");

        ALPAKA_FN_HOST_ACC constexpr RayWalkBehaviour(T_Terms... terms) : T_Terms{terms}...
        {
        }

        ALPAKA_FN_ACC bool onCell(
            auto const& acc,
            hase::core::DeviceMeshView const& mesh,
            auto& rayState,
            unsigned const cell,
            Tet4FaceIntersection const intersection)
        {
            return (invokeCell(static_cast<T_Terms&>(*this), acc, mesh, rayState, cell, intersection) && ...);
        }

        ALPAKA_FN_ACC BoundaryResult onBoundary(
            auto const& acc,
            hase::core::DeviceMeshView const& mesh,
            auto& rayState,
            unsigned const cell,
            unsigned const localFace)
        {
            bool const continueTraversal
                = ((invokeBoundary(static_cast<T_Terms&>(*this), acc, mesh, rayState, cell, localFace)
                    == BoundaryResult::continueTraversal)
                   && ...);
            return continueTraversal ? BoundaryResult::continueTraversal : BoundaryResult::stop;
        }

    private:
        template<typename T_Term, typename T_Acc, typename T_RayState>
        requires HasOnCell<T_Term, T_Acc, T_RayState>
        ALPAKA_FN_ACC static bool invokeCell(
            T_Term& term,
            T_Acc const& acc,
            hase::core::DeviceMeshView const& mesh,
            T_RayState& rayState,
            unsigned const cell,
            Tet4FaceIntersection const intersection)
        {
            return term(acc, mesh, rayState, cell, intersection);
        }

        template<typename T_Term, typename T_Acc, typename T_RayState>
        requires(!HasOnCell<T_Term, T_Acc, T_RayState>)
        ALPAKA_FN_ACC static bool invokeCell(
            T_Term&,
            T_Acc const&,
            hase::core::DeviceMeshView const&,
            T_RayState&,
            unsigned,
            Tet4FaceIntersection)
        {
            return true;
        }

        template<typename T_Term, typename T_Acc, typename T_RayState>
        requires HasOnBoundary<T_Term, T_Acc, T_RayState>
        ALPAKA_FN_ACC static BoundaryResult invokeBoundary(
            T_Term& term,
            T_Acc const& acc,
            hase::core::DeviceMeshView const& mesh,
            T_RayState& rayState,
            unsigned const cell,
            unsigned const localFace)
        {
            return term(acc, mesh, rayState, cell, localFace);
        }

        template<typename T_Term, typename T_Acc, typename T_RayState>
        requires(!HasOnBoundary<T_Term, T_Acc, T_RayState>)
        ALPAKA_FN_ACC static BoundaryResult invokeBoundary(
            T_Term&,
            T_Acc const&,
            hase::core::DeviceMeshView const&,
            T_RayState&,
            unsigned,
            unsigned)
        {
            return BoundaryResult::continueTraversal;
        }
    };

    template<typename... T_Terms>
    RayWalkBehaviour(T_Terms...) -> RayWalkBehaviour<T_Terms...>;

    /**
     * Walk one ray through the Tet4 mesh.
     *
     * Cell contribution and boundary handling are compile-time policies. The
     * only branches left here are geometric state transitions that every ray
     * tracer must perform.
     */
    enum class WalkResult : std::uint8_t
    {
        reachedBoundary,
        terminatedInCell,
        failed
    };

    template<State T_Ray, typename T_Behaviour>
    [[nodiscard]] ALPAKA_FN_ACC WalkResult walk(
        auto const& acc,
        hase::core::DeviceMeshView const& mesh,
        T_Ray& rayState,
        T_Behaviour behaviour)
    {
        constexpr unsigned maxTraversalSteps = 10000u;
        for(unsigned step = 0u; step < maxTraversalSteps; ++step)
        {
            unsigned const currentCell = rayState.cell;
            assert(currentCell < mesh.numberOfCells);
            auto const intersection = nextFaceIntersection(
                mesh,
                currentCell,
                rayState.position,
                rayState.direction,
                rayState.forbiddenFace);
            if(intersection.localFace < 0)
            {
                int const recoveryFace = isNearTet4Face(mesh, currentCell, rayState.position)
                                             ? immediateExitFace(
                                                   mesh,
                                                   currentCell,
                                                   rayState.position,
                                                   rayState.direction,
                                                   rayState.forbiddenFace)
                                             : -1;
                if(recoveryFace < 0)
                {
                    return WalkResult::failed;
                }
                auto const transition
                    = recoverFaceTransition(mesh, currentCell, recoveryFace, rayState.position, rayState.direction);
                if(transition.status == Tet4TransitionStatus::failed)
                {
                    rayState.cell = transition.cell;
                    return WalkResult::failed;
                }
                if(transition.status == Tet4TransitionStatus::reachedBoundary)
                {
                    rayState.cell = transition.cell;
                    auto const boundaryResult = behaviour.onBoundary(
                        acc,
                        mesh,
                        rayState,
                        transition.cell,
                        static_cast<unsigned>(transition.boundaryFace));
                    if(boundaryResult == BoundaryResult::continueTraversal)
                        continue;
                    return WalkResult::reachedBoundary;
                }
                rayState.cell = transition.cell;
                rayState.forbiddenFace = transition.forbiddenFace;
                continue;
            }

            if(!behaviour.onCell(acc, mesh, rayState, currentCell, intersection))
            {
                return WalkResult::terminatedInCell;
            }
            rayState.position = advance(rayState.position, rayState.direction, intersection.length);
            auto const transition
                = transitionAcrossIntersection(mesh, currentCell, intersection, rayState.position, rayState.direction);
            if(transition.status == Tet4TransitionStatus::failed)
            {
                rayState.cell = transition.cell;
                return WalkResult::failed;
            }
            if(transition.status == Tet4TransitionStatus::reachedBoundary)
            {
                rayState.cell = transition.cell;
                auto const boundaryResult = behaviour.onBoundary(
                    acc,
                    mesh,
                    rayState,
                    transition.cell,
                    static_cast<unsigned>(transition.boundaryFace));
                if(boundaryResult == BoundaryResult::continueTraversal)
                    continue;
                return WalkResult::reachedBoundary;
            }
            rayState.cell = transition.cell;
            rayState.forbiddenFace = transition.forbiddenFace;
        }
        return WalkResult::failed;
    }
} // namespace hase::kernels::forward::ray
