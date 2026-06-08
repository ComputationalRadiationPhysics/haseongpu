/* Copyright 2026 Tim Hanel
 */
#pragma once
#include <cstdint>
#include <random>

namespace hase::random
{
    namespace internal
    {
        inline constexpr std::uint64_t boostHashMix(std::uint64_t x)
        {
            std::uint64_t const m = 0xe9'846a'fb1a'615dull;
            x ^= x >> 32;
            x *= m;
            x ^= x >> 32;
            x *= m;
            x ^= x >> 28;
            return x;
        }

        inline constexpr std::uint32_t mixSeed(std::uint32_t seed, std::uint32_t value)
        {
            // boost::hash_combine
            return boostHashMix(seed + 0x9e37'79b9 + value);
        }
    } // namespace internal

    struct SeedGenerator
    {
        // Default non-deterministic seed generated.
        explicit SeedGenerator() : fixedSeed(std::random_device{}())
        {
        }

        SeedGenerator(SeedGenerator const&) = delete;
        SeedGenerator& operator=(SeedGenerator const&) = delete;

        SeedGenerator(SeedGenerator&&) = delete;
        SeedGenerator& operator=(SeedGenerator&&) = delete;

        // allows setting a custom fixed seed from the outside interface
        void updateSeed(unsigned const seed)
        {
            fixedSeed = seed;
        }

        [[nodiscard]] unsigned getSeed() const
        {
            return fixedSeed;
        }

        static SeedGenerator& get()
        {
            static SeedGenerator provider{};
            return provider;
        }

    private:
        unsigned fixedSeed;
    };

    inline constexpr std::uint32_t seedForWorker(std::uint32_t base, std::uint32_t rank, std::uint32_t deviceIndex)
    {
        return internal::mixSeed(internal::mixSeed(base, rank), deviceIndex);
    }

} // namespace hase::random
