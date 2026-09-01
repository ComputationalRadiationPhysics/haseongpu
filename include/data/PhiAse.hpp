#pragma once

#include <transport/TransportReader.hpp>

#include <cstdint>
#include <optional>
#include <string>

namespace hase::data
{
    /** @brief User-facing ASE trace controls independent of resident device data. */
    class PhiAse
    {
    public:
        struct FieldName
        {
            static constexpr char const* propagationMode = "propagationMode";
            static constexpr char const* minRays = "minRays";
            static constexpr char const* maxRays = "maxRays";
            static constexpr char const* forwardRayCount = "forwardRayCount";
            static constexpr char const* relativeStandardErrorThreshold = "relativeStandardErrorThreshold";
            static constexpr char const* trackRayVisits = "trackRayVisits";
            static constexpr char const* repetitions = "repetitions";
            static constexpr char const* adaptiveSteps = "adaptiveSteps";
            static constexpr char const* useReflections = "useReflections";
            static constexpr char const* reflectionMode = "reflectionMode";
            static constexpr char const* surfaceReservoirSize = "surfaceReservoirSize";
            static constexpr char const* srmPositionMode = "srmPositionMode";
            static constexpr char const* reflectionMaxIterations = "reflectionMaxIterations";
            static constexpr char const* reflectionTolerance = "reflectionTolerance";
            static constexpr char const* monochromatic = "monochromatic";
            static constexpr char const* backend = "backend";
            static constexpr char const* parallelMode = "parallelMode";
            static constexpr char const* numDevices = "numDevices";
            static constexpr char const* writeVtk = "writeVtk";
            static constexpr char const* devices = "devices";
            static constexpr char const* minSampleRange = "minSampleRange";
            static constexpr char const* maxSampleRange = "maxSampleRange";
            static constexpr char const* rngSeed = "rngSeed";
            static constexpr char const* aseSteps = "aseSteps";
        };

        std::string propagationMode;
        std::uint64_t minRays{};
        std::uint64_t maxRays{};
        std::optional<std::uint64_t> forwardRayCount;
        double relativeStandardErrorThreshold{};
        bool trackRayVisits{};
        std::uint64_t repetitions{};
        std::uint64_t adaptiveSteps{};
        bool useReflections{};
        std::string reflectionMode;
        std::uint64_t surfaceReservoirSize{};
        std::string srmPositionMode;
        std::uint64_t reflectionMaxIterations{};
        double reflectionTolerance{};
        bool monochromatic{};
        std::optional<std::string> backend;
        std::string parallelMode;
        std::uint64_t numDevices{};
        bool writeVtk{};
        std::optional<transport::Array<std::uint64_t>> devices;
        std::optional<std::uint64_t> minSampleRange;
        std::optional<std::uint64_t> maxSampleRange;
        std::optional<std::uint64_t> rngSeed;
        std::optional<std::uint64_t> aseSteps;

        /**
         * @brief Read user-facing ASE controls from one graph node.
         * @param reader Typed reader for the active transport iteration.
         * @param prefix Path of the PhiASE node.
         * @return Transported ASE controls without execution preparation.
         */
        static PhiAse fromTransport(transport::TransportReader const& reader, transport::TransportPath const& prefix);
    };
} // namespace hase::data
