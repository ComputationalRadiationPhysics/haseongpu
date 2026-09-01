/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <alpaka/alpaka.hpp>

#include <concepts>
#include <type_traits>

namespace hase::kernels::forward::tracePolicy
{
    struct Policy
    {
    };

    namespace category
    {
        struct Source : Policy
        {
            template<typename T_Policy>
            constexpr bool operator()(T_Policy) const
            {
                return std::derived_from<std::remove_cvref_t<T_Policy>, Source>;
            }
        };

        struct Cell : Policy
        {
            template<typename T_Policy>
            constexpr bool operator()(T_Policy) const
            {
                return std::derived_from<std::remove_cvref_t<T_Policy>, Cell>;
            }
        };

        struct Boundary : Policy
        {
            template<typename T_Policy>
            constexpr bool operator()(T_Policy) const
            {
                return std::derived_from<std::remove_cvref_t<T_Policy>, Boundary>;
            }
        };

        struct Position : Policy
        {
            template<typename T_Policy>
            constexpr bool operator()(T_Policy) const
            {
                return std::derived_from<std::remove_cvref_t<T_Policy>, Position>;
            }
        };

        struct Diagnostics : Policy
        {
            template<typename T_Policy>
            constexpr bool operator()(T_Policy) const
            {
                return std::derived_from<std::remove_cvref_t<T_Policy>, Diagnostics>;
            }
        };
    } // namespace category

    namespace source
    {
        struct Volume : category::Source
        {
        };

        struct ReflectionCandidates : category::Source
        {
        };

        struct SurfaceReservoir : category::Source
        {
        };

        inline constexpr Volume volume;
        inline constexpr ReflectionCandidates reflectionCandidates;
        inline constexpr SurfaceReservoir surfaceReservoir;
    } // namespace source

    namespace cell
    {
        struct ForwardAse : category::Cell
        {
        };

        inline constexpr ForwardAse forwardAse;
    } // namespace cell

    namespace boundary
    {
        struct Escape : category::Boundary
        {
        };

        struct ReflectionCandidates : category::Boundary
        {
        };

        struct SurfaceReservoir : category::Boundary
        {
        };

        inline constexpr Escape escape;
        inline constexpr ReflectionCandidates reflectionCandidates;
        inline constexpr SurfaceReservoir surfaceReservoir;
    } // namespace boundary

    namespace position
    {
        struct None : category::Position
        {
        };

        struct Exact : category::Position
        {
        };

        struct Centroid : category::Position
        {
        };

        inline constexpr None none;
        inline constexpr Exact exact;
        inline constexpr Centroid centroid;
    } // namespace position

    namespace diagnostics
    {
        struct None : category::Diagnostics
        {
        };

        struct CellRayVisits : category::Diagnostics
        {
        };

        inline constexpr None none;
        inline constexpr CellRayVisits cellRayVisits;
    } // namespace diagnostics
} // namespace hase::kernels::forward::tracePolicy

namespace hase::kernels::forward::concepts
{
    template<typename T_Policy>
    concept TracePolicy = std::derived_from<std::remove_cvref_t<T_Policy>, tracePolicy::Policy>
                          && std::default_initializable<std::remove_cvref_t<T_Policy>>;
} // namespace hase::kernels::forward::concepts

namespace alpaka::trait
{
    template<hase::kernels::forward::concepts::TracePolicy T_Policy>
    struct IsPolicy<T_Policy> : std::true_type
    {
    };
} // namespace alpaka::trait

namespace hase::kernels::forward
{
    /** @brief Unordered, stateless compile-time configuration for one forward trace kernel. */
    template<concepts::TracePolicy... T_Policies>
    struct TracePolicyList : alpaka::PolicyList<T_Policies...>
    {
        using Base = alpaka::PolicyList<T_Policies...>;

        constexpr TracePolicyList(T_Policies... policies) : Base{policies...}
        {
        }

        static consteval auto getSource()
        {
            return Base::search(tracePolicy::category::Source{}, tracePolicy::source::volume);
        }

        static consteval auto getCell()
        {
            return Base::search(tracePolicy::category::Cell{}, tracePolicy::cell::forwardAse);
        }

        static consteval auto getBoundary()
        {
            return Base::search(tracePolicy::category::Boundary{}, tracePolicy::boundary::escape);
        }

        static consteval auto getPosition()
        {
            return Base::search(tracePolicy::category::Position{}, tracePolicy::position::none);
        }

        static consteval auto getDiagnostics()
        {
            return Base::search(tracePolicy::category::Diagnostics{}, tracePolicy::diagnostics::none);
        }

        using Base::hasPolicy;
    };

    template<typename... T_Policies>
    TracePolicyList(T_Policies...) -> TracePolicyList<T_Policies...>;
} // namespace hase::kernels::forward

namespace hase::kernels::forward::concepts
{
    template<typename T_Policies>
    concept TracePolicyList = alpaka::concepts::SpecializationOf<T_Policies, hase::kernels::forward::TracePolicyList>;
} // namespace hase::kernels::forward::concepts
