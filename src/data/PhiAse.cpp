#include <data/PhiAse.hpp>

namespace hase::data
{
    PhiAse PhiAse::fromTransport(transport::TransportReader const& reader, transport::TransportPath const& prefix)
    {
        reader.prefetch(prefix);
        PhiAse result;
        reader.assign(result.propagationMode, prefix, FieldName::propagationMode);
        reader.assign(result.minRays, prefix, FieldName::minRays);
        reader.assign(result.maxRays, prefix, FieldName::maxRays);
        reader.assign(result.forwardRayCount, prefix, FieldName::forwardRayCount);
        reader.assign(result.relativeStandardErrorThreshold, prefix, FieldName::relativeStandardErrorThreshold);
        if(reader.contains(prefix, FieldName::trackRayVisits))
            reader.assign(result.trackRayVisits, prefix, FieldName::trackRayVisits);
        reader.assign(result.repetitions, prefix, FieldName::repetitions);
        reader.assign(result.adaptiveSteps, prefix, FieldName::adaptiveSteps);
        reader.assign(result.useReflections, prefix, FieldName::useReflections);
        reader.assign(result.reflectionMode, prefix, FieldName::reflectionMode);
        reader.assign(result.surfaceReservoirSize, prefix, FieldName::surfaceReservoirSize);
        reader.assign(result.srmPositionMode, prefix, FieldName::srmPositionMode);
        reader.assign(result.reflectionMaxIterations, prefix, FieldName::reflectionMaxIterations);
        reader.assign(result.reflectionTolerance, prefix, FieldName::reflectionTolerance);
        reader.assign(result.monochromatic, prefix, FieldName::monochromatic);
        reader.assign(result.backend, prefix, FieldName::backend);
        reader.assign(result.parallelMode, prefix, FieldName::parallelMode);
        reader.assign(result.numDevices, prefix, FieldName::numDevices);
        reader.assign(result.writeVtk, prefix, FieldName::writeVtk);
        reader.assign(result.devices, prefix, FieldName::devices);
        reader.assign(result.minSampleRange, prefix, FieldName::minSampleRange);
        reader.assign(result.maxSampleRange, prefix, FieldName::maxSampleRange);
        reader.assign(result.rngSeed, prefix, FieldName::rngSeed);
        reader.assign(result.aseSteps, prefix, FieldName::aseSteps);
        return result;
    }
} // namespace hase::data
