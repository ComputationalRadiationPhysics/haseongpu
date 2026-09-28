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
#include <cstdint>
#include <span>
#include <type_traits>

namespace hase::kernels
{
    /** @brief Element-wise `base + scale * slope` update. */
    struct AddScaled
    {
        double scale = 1.0;

        ALPAKA_FN_ACC constexpr auto operator()(
            alpaka::concepts::Simd auto const& base,
            alpaka::concepts::Simd auto const& slope) const
        {
            return base + scale * slope;
        }
    };

    /** @brief Element-wise second-order Heun stage combination. */
    struct CombineHeun
    {
        double timeStep = 0.0;

        ALPAKA_FN_ACC constexpr auto operator()(
            alpaka::concepts::Simd auto const& base,
            alpaka::concepts::Simd auto const& first,
            alpaka::concepts::Simd auto const& second) const
        {
            return base + 0.5 * timeStep * (first + second);
        }
    };

    /** @brief Element-wise classical fourth-order Runge-Kutta combination. */
    struct CombineRungeKutta4
    {
        double timeStep = 0.0;

        ALPAKA_FN_ACC constexpr auto operator()(
            alpaka::concepts::Simd auto const& base,
            alpaka::concepts::Simd auto const& k1,
            alpaka::concepts::Simd auto const& k2,
            alpaka::concepts::Simd auto const& k3,
            alpaka::concepts::Simd auto const& k4) const
        {
            return base + (timeStep / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }
    };

    /** @brief Compose a stage derivative while reusing frozen pump and ASE rates. */
    struct ComposeFrozenSourcesDerivative
    {
        std::span<std::uint8_t const> materialActive;
        std::span<double const> materialFluorescenceLifetimes;

        ALPAKA_FN_ACC auto operator()(
            alpaka::concepts::Simd auto const& betaVolume,
            alpaka::concepts::Simd auto const& dndtPump,
            alpaka::concepts::Simd auto const& dndtAse,
            alpaka::concepts::Simd auto const& materialIds) const
        {
            using T_Result = std::remove_cvref_t<decltype(betaVolume)>;
            static_assert(std::same_as<alpaka::GetValueType_t<T_Result>, double>);
            static_assert(std::same_as<alpaka::GetValueType_t<std::remove_cvref_t<decltype(materialIds)>>, unsigned>);
            auto const active = alpaka::SimdMask<double, T_Result::width()>{
                [&](auto const lane)
                {
                    auto const material = materialIds[decltype(lane)::value];
                    return materialActive[material] != 0u;
                }};
            auto const tau
                = T_Result{[&](auto const lane)
                           {
                               auto const material = materialIds[decltype(lane)::value];
                               return active[decltype(lane)::value] ? materialFluorescenceLifetimes[material] : 1.0;
                           }};
            auto const derivative = dndtPump - dndtAse - betaVolume / tau;
            auto result = T_Result{[](auto) { return 0.0; }};
            alpaka::where(active, result) = derivative;
            return result;
        }
    };

    /**
     * @brief Exact fluorescence-decay update with frozen pump and ASE sources.
     *
     * Material activity and lifetime are gathered for each SIMD lane through
     * its cell material id. Inactive lanes are masked to exactly zero.
     */
    struct ExponentialEulerUpdate
    {
        double timeStep = 0.0;
        std::span<std::uint8_t const> materialActive;
        std::span<double const> materialFluorescenceLifetimes;

        /**
         * @param betaVolume Cell excitation fractions.
         * @param dndtPump Frozen pump source rates.
         * @param dndtAse Frozen ASE depletion rates.
         * @param materialIds Cell material ids used for property gathering.
         * @return Updated excitation fractions, with inactive lanes set to zero.
         */
        ALPAKA_FN_ACC auto operator()(
            alpaka::concepts::Simd auto const& betaVolume,
            alpaka::concepts::Simd auto const& dndtPump,
            alpaka::concepts::Simd auto const& dndtAse,
            alpaka::concepts::Simd auto const& materialIds) const
        {
            using T_Result = std::remove_cvref_t<decltype(betaVolume)>;
            static_assert(std::same_as<alpaka::GetValueType_t<T_Result>, double>);
            static_assert(std::same_as<alpaka::GetValueType_t<std::remove_cvref_t<decltype(materialIds)>>, unsigned>);
            auto const active = alpaka::SimdMask<double, T_Result::width()>{
                [&](auto const lane)
                {
                    auto const material = materialIds[decltype(lane)::value];
                    return materialActive[material] != 0u;
                }};
            auto const tau
                = T_Result{[&](auto const lane)
                           {
                               auto const material = materialIds[decltype(lane)::value];
                               return active[decltype(lane)::value] ? materialFluorescenceLifetimes[material] : 1.0;
                           }};
            auto const decay = alpaka::math::exp(-timeStep / tau);
            auto const source = dndtPump - dndtAse;
            auto const updated = tau * source * (1.0 - decay) + betaVolume * decay;
            auto result = T_Result{[](auto) { return 0.0; }};
            alpaka::where(active, result) = updated;
            return result;
        }
    };

    /** @brief Clamp excitation fractions and map NaN values into `[0, 1]`. */
    struct ClipBeta
    {
        ALPAKA_FN_ACC constexpr auto operator()(alpaka::concepts::Simd auto const& betaVolume) const
        {
            auto result = betaVolume;
            alpaka::where((result <= 0.0) || (result != result), result) = 0.0;
            alpaka::where(result > 1.0, result) = 1.0;
            return result;
        }
    };

} // namespace hase::kernels
