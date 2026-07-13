/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <alpakaUtils/DevBundle.hpp>
#include <alpakaUtils/memory.hpp>
#include <concepts/concepts.hpp>
#include <core/forwardSrm.hpp>
#include <core/mesh.hpp>

#include <cstddef>
#include <utility>
#include <vector>

namespace hase::core
{
    namespace detail
    {
        template<typename T_Buffer>
        std::vector<typename T_Buffer::value_type> copyBufferToVector(
            concepts::Queue auto const& queue,
            T_Buffer const& buffer)
        {
            auto hostBuffer = alpaka::onHost::allocHostLike(buffer);
            alpaka::onHost::memcpy(queue, hostBuffer, buffer);
            alpaka::onHost::wait(queue);
            auto const* data = alpaka::onHost::data(hostBuffer);
            return {data, data + buffer.getExtents().product()};
        }

    } // namespace detail

    template<alpaka::onHost::concepts::Device T_Device, alpaka::concepts::Executor T_Exec>
    class ForwardSrmDeviceState
    {
        using T_Queue = ALPAKA_TYPEOF(std::declval<T_Device>().makeQueue(alpaka::queueKind::blocking));
        using T_DoubleBuffer = ALPAKA_TYPEOF(alpaka::onHost::alloc<double>(std::declval<T_Device&>(), std::size_t{1}));
        using T_UnsignedBuffer
            = ALPAKA_TYPEOF(alpaka::onHost::alloc<unsigned>(std::declval<T_Device&>(), std::size_t{1}));
        using T_CharBuffer = ALPAKA_TYPEOF(alpaka::onHost::alloc<char>(std::declval<T_Device&>(), std::size_t{1}));
        using T_Accumulation = hase::kernels::forward::ForwardAccumulationSpans<
            T_DoubleBuffer,
            T_DoubleBuffer,
            T_UnsignedBuffer,
            T_UnsignedBuffer>;
        using T_Spectrum = hase::kernels::forward::ForwardSpectrumSpans<T_DoubleBuffer, T_DoubleBuffer>;
        using T_Reservoir = hase::kernels::forward::SurfaceReservoirSpans<
            T_UnsignedBuffer,
            T_DoubleBuffer,
            T_DoubleBuffer,
            T_DoubleBuffer,
            T_DoubleBuffer,
            T_UnsignedBuffer,
            T_DoubleBuffer>;
        using T_SamplingCdf = hase::kernels::forward::SurfaceReservoirSamplingCdfSpans<
            T_DoubleBuffer,
            T_DoubleBuffer,
            T_UnsignedBuffer>;

    public:
        ForwardSrmDeviceState(
            T_Device const& device,
            T_Exec const& executor,
            DeviceMeshContainer<T_Device> const& meshContainer,
            ExperimentParameters const& experiment,
            double const betaVolumeTotal,
            unsigned const rayCount,
            unsigned const rngSeed)
            : m_devBundle(device, executor)
            , m_queue(m_devBundle.device.makeQueue(alpaka::queueKind::blocking))
            , m_mesh(meshContainer.toView())
            , m_rayCount(rayCount)
            , m_rngSeed(rngSeed)
            , m_betaVolumeTotal(betaVolumeTotal)
            , m_slotsPerFace(experiment.surfaceReservoirSize)
            , m_phi(alpaka::onHost::alloc<double>(m_devBundle.device, static_cast<std::size_t>(m_mesh.numberOfCells)))
            , m_phiSquare(
                  alpaka::onHost::alloc<double>(m_devBundle.device, static_cast<std::size_t>(m_mesh.numberOfCells)))
            , m_volumeRayVisits(
                  alpaka::onHost::alloc<unsigned>(m_devBundle.device, static_cast<std::size_t>(m_mesh.numberOfCells)))
            , m_droppedRays(
                  alpaka::onHost::alloc<unsigned>(m_devBundle.device, static_cast<std::size_t>(m_mesh.numberOfCells)))
            , m_sigmaA(hase::alpakaUtils::toDevice(m_queue, experiment.sigmaA))
            , m_sigmaE(hase::alpakaUtils::toDevice(m_queue, experiment.sigmaE))
            , m_accumulation{m_phi, m_phiSquare, m_volumeRayVisits, m_droppedRays}
            , m_spectrum{m_sigmaA, m_sigmaE, static_cast<unsigned>(experiment.sigmaA.size())}
            , m_countsA(alpaka::onHost::alloc<unsigned>(m_devBundle.device, static_cast<std::size_t>(faceCount())))
            , m_countsB(alpaka::onHost::alloc<unsigned>(m_devBundle.device, static_cast<std::size_t>(faceCount())))
            , m_dirXA(alpaka::onHost::alloc<double>(m_devBundle.device, static_cast<std::size_t>(reservoirSlots())))
            , m_dirXB(alpaka::onHost::alloc<double>(m_devBundle.device, static_cast<std::size_t>(reservoirSlots())))
            , m_dirYA(alpaka::onHost::alloc<double>(m_devBundle.device, static_cast<std::size_t>(reservoirSlots())))
            , m_dirYB(alpaka::onHost::alloc<double>(m_devBundle.device, static_cast<std::size_t>(reservoirSlots())))
            , m_dirZA(alpaka::onHost::alloc<double>(m_devBundle.device, static_cast<std::size_t>(reservoirSlots())))
            , m_dirZB(alpaka::onHost::alloc<double>(m_devBundle.device, static_cast<std::size_t>(reservoirSlots())))
            , m_weightsA(alpaka::onHost::alloc<double>(m_devBundle.device, static_cast<std::size_t>(reservoirSlots())))
            , m_weightsB(alpaka::onHost::alloc<double>(m_devBundle.device, static_cast<std::size_t>(reservoirSlots())))
            , m_sigmaIndicesA(
                  alpaka::onHost::alloc<unsigned>(m_devBundle.device, static_cast<std::size_t>(reservoirSlots())))
            , m_sigmaIndicesB(
                  alpaka::onHost::alloc<unsigned>(m_devBundle.device, static_cast<std::size_t>(reservoirSlots())))
            , m_faceWeightsA(
                  alpaka::onHost::alloc<double>(m_devBundle.device, static_cast<std::size_t>(faceCount())))
            , m_faceWeightsB(
                  alpaka::onHost::alloc<double>(m_devBundle.device, static_cast<std::size_t>(faceCount())))
            , m_samplingCdf(alpaka::onHost::alloc<double>(m_devBundle.device, static_cast<std::size_t>(faceCount())))
            , m_samplingTotalWeight(alpaka::onHost::alloc<double>(m_devBundle.device, std::size_t{1}))
            , m_systematicOffset(alpaka::onHost::alloc<double>(m_devBundle.device, std::size_t{1}))
            , m_stratifiedRayCounts(
                  alpaka::onHost::alloc<unsigned>(m_devBundle.device, static_cast<std::size_t>(faceCount())))
            , m_stratifiedRayOffsets(
                  alpaka::onHost::alloc<unsigned>(m_devBundle.device, static_cast<std::size_t>(faceCount())))
            , m_stratifiedRayFaces(
                  alpaka::onHost::alloc<unsigned>(m_devBundle.device, static_cast<std::size_t>(m_rayCount)))
            , m_samplingCdfScanBuffer(alpaka::onHost::alloc<char>(
                  m_devBundle.device,
                  alpaka::onHost::getScanBufferSize<double>(alpaka::Vec{static_cast<std::size_t>(faceCount())})))
            , m_stratifiedCountScanBuffer(alpaka::onHost::alloc<char>(
                  m_devBundle.device,
                  alpaka::onHost::getScanBufferSize<unsigned>(alpaka::Vec{static_cast<std::size_t>(faceCount())})))
            , m_reservoirA{
                  m_countsA,
                  m_dirXA,
                  m_dirYA,
                  m_dirZA,
                  m_weightsA,
                  m_sigmaIndicesA,
                  m_faceWeightsA,
                  m_slotsPerFace}
            , m_reservoirB{
                  m_countsB,
                  m_dirXB,
                  m_dirYB,
                  m_dirZB,
                  m_weightsB,
                  m_sigmaIndicesB,
                  m_faceWeightsB,
                  m_slotsPerFace}
            , m_samplingCdfSpans{
                  m_samplingCdf,
                  m_samplingTotalWeight,
                  m_stratifiedRayFaces,
                  faceCount() <= m_rayCount}
        {
            alpaka::onHost::fill(
                m_queue,
                m_phi,
                0.0,
                alpaka::Vec{static_cast<std::size_t>(m_mesh.numberOfCells)});
            alpaka::onHost::fill(
                m_queue,
                m_phiSquare,
                0.0,
                alpaka::Vec{static_cast<std::size_t>(m_mesh.numberOfCells)});
            alpaka::onHost::fill(
                m_queue,
                m_volumeRayVisits,
                0u,
                alpaka::Vec{static_cast<std::size_t>(m_mesh.numberOfCells)});
            alpaka::onHost::fill(
                m_queue,
                m_droppedRays,
                0u,
                alpaka::Vec{static_cast<std::size_t>(m_mesh.numberOfCells)});
            clearReservoir(m_reservoirA);
            clearReservoir(m_reservoirB);
            alpaka::onHost::wait(m_queue);
        }

        void traceDirect()
        {
            auto const frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                m_devBundle.device,
                m_devBundle.executor,
                alpaka::Vec{m_rayCount});
            m_queue.enqueue(
                frameSpec,
                alpaka::KernelBundle{
                    hase::kernels::forward::AccumulateForwardPhiAseReservoir{},
                    m_mesh,
                    m_rayCount,
                    m_betaVolumeTotal,
                    m_accumulation,
                    m_reservoirA,
                    m_spectrum,
                    m_rngSeed});
            alpaka::onHost::wait(m_queue);
            m_inputReservoirA = true;
        }

        [[nodiscard]] double prepareDirectSamplingCdf()
        {
            return prepareSamplingCdf(m_reservoirA, m_rngSeed);
        }

        [[nodiscard]] unsigned rayCount() const
        {
            return m_rayCount;
        }

        [[nodiscard]] double prepareReflectedSamplingCdf(unsigned const pass)
        {
            return prepareSamplingCdf(currentInputReservoir(), m_rngSeed + pass * m_rayCount);
        }

        void traceReflected(double const sourceWeight, unsigned const pass)
        {
            T_Reservoir& inputReservoir = currentInputReservoir();
            T_Reservoir& outputReservoir = m_inputReservoirA ? m_reservoirB : m_reservoirA;
            clearReservoir(outputReservoir);
            alpaka::onHost::wait(m_queue);
            if(sourceWeight <= 0.0)
            {
                m_inputReservoirA = !m_inputReservoirA;
                return;
            }
            auto const frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                m_devBundle.device,
                m_devBundle.executor,
                alpaka::Vec{m_rayCount});
            m_queue.enqueue(
                frameSpec,
                alpaka::KernelBundle{
                    hase::kernels::forward::AccumulateReflectedForwardPhiAse{},
                    m_mesh,
                    m_rayCount,
                    sourceWeight,
                    m_accumulation,
                    inputReservoir,
                    m_samplingCdfSpans,
                    outputReservoir,
                    m_spectrum,
                    m_rngSeed + pass * m_rayCount});
            alpaka::onHost::wait(m_queue);
            m_inputReservoirA = !m_inputReservoirA;
        }

        void downloadAccumulation(ForwardPhiAseRawResult& result)
        {
            result.scoreSum = detail::copyBufferToVector(m_queue, m_phi);
            result.scoreSquareSum = detail::copyBufferToVector(m_queue, m_phiSquare);
            result.totalRays = detail::copyBufferToVector(m_queue, m_volumeRayVisits);
            result.droppedRays = detail::copyBufferToVector(m_queue, m_droppedRays);
            result.rayCount = m_rayCount;
        }

    private:
        [[nodiscard]] unsigned faceCount() const
        {
            return m_mesh.numberOfCells * m_mesh.numberOfFacesPerCell;
        }

        [[nodiscard]] unsigned reservoirSlots() const
        {
            return faceCount() * m_slotsPerFace;
        }

        [[nodiscard]] T_Reservoir& currentInputReservoir()
        {
            return m_inputReservoirA ? m_reservoirA : m_reservoirB;
        }

        [[nodiscard]] double prepareSamplingCdf(T_Reservoir const& reservoir, unsigned const seed)
        {
            auto const faceFrameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                m_devBundle.device,
                m_devBundle.executor,
                alpaka::Vec{faceCount()});
            auto const scalarFrameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                m_devBundle.device,
                m_devBundle.executor,
                alpaka::Vec{1u});
            alpaka::onHost::inclusiveScan(
                m_queue,
                m_devBundle.executor,
                m_samplingCdfScanBuffer,
                m_samplingCdf,
                reservoir.faceWeights);
            m_queue.enqueue(
                scalarFrameSpec,
                alpaka::KernelBundle{
                    hase::kernels::forward::CaptureSurfaceReservoirSamplingTotalWeight{},
                    faceCount(),
                    m_samplingCdfSpans});
            m_queue.enqueue(
                faceFrameSpec,
                alpaka::KernelBundle{
                    hase::kernels::forward::NormalizeSurfaceReservoirSamplingCdf{},
                    faceCount(),
                    m_samplingCdfSpans});
            if(m_samplingCdfSpans.useFaceStratification)
            {
                m_queue.enqueue(
                    scalarFrameSpec,
                    alpaka::KernelBundle{
                        hase::kernels::forward::GenerateSurfaceReservoirSystematicOffset{},
                        m_systematicOffset,
                        seed});
                m_queue.enqueue(
                    faceFrameSpec,
                    alpaka::KernelBundle{
                        hase::kernels::forward::AssignSurfaceReservoirStratifiedRayCounts{},
                        faceCount(),
                        m_rayCount,
                        m_samplingCdfSpans,
                        m_systematicOffset,
                        m_stratifiedRayCounts});
                alpaka::onHost::exclusiveScan(
                    m_queue,
                    m_devBundle.executor,
                    m_stratifiedCountScanBuffer,
                    m_stratifiedRayOffsets,
                    m_stratifiedRayCounts);
                m_queue.enqueue(
                    faceFrameSpec,
                    alpaka::KernelBundle{
                        hase::kernels::forward::ScatterSurfaceReservoirStratifiedRayFaces{},
                        faceCount(),
                        m_stratifiedRayCounts,
                        m_stratifiedRayOffsets,
                        m_stratifiedRayFaces});
            }
            alpaka::onHost::wait(m_queue);

            auto hostTotalWeight = alpaka::onHost::allocHostLike(m_samplingTotalWeight);
            alpaka::onHost::memcpy(m_queue, hostTotalWeight, m_samplingTotalWeight);
            alpaka::onHost::wait(m_queue);
            return alpaka::onHost::data(hostTotalWeight)[0u];
        }

        void clearReservoir(T_Reservoir& reservoir)
        {
            alpaka::onHost::fill(
                m_queue,
                reservoir.counts,
                0u,
                alpaka::Vec{static_cast<std::size_t>(faceCount())});
            alpaka::onHost::fill(
                m_queue,
                reservoir.faceWeights,
                0.0,
                alpaka::Vec{static_cast<std::size_t>(faceCount())});
        }

        hase::alpakaUtils::DevBundle<T_Device, T_Exec> m_devBundle;
        T_Queue m_queue;
        DeviceMeshView m_mesh;
        unsigned m_rayCount;
        unsigned m_rngSeed;
        double m_betaVolumeTotal;
        unsigned m_slotsPerFace;
        T_DoubleBuffer m_phi;
        T_DoubleBuffer m_phiSquare;
        T_UnsignedBuffer m_volumeRayVisits;
        T_UnsignedBuffer m_droppedRays;
        T_DoubleBuffer m_sigmaA;
        T_DoubleBuffer m_sigmaE;
        T_Accumulation m_accumulation;
        T_Spectrum m_spectrum;
        T_UnsignedBuffer m_countsA;
        T_UnsignedBuffer m_countsB;
        T_DoubleBuffer m_dirXA;
        T_DoubleBuffer m_dirXB;
        T_DoubleBuffer m_dirYA;
        T_DoubleBuffer m_dirYB;
        T_DoubleBuffer m_dirZA;
        T_DoubleBuffer m_dirZB;
        T_DoubleBuffer m_weightsA;
        T_DoubleBuffer m_weightsB;
        T_UnsignedBuffer m_sigmaIndicesA;
        T_UnsignedBuffer m_sigmaIndicesB;
        T_DoubleBuffer m_faceWeightsA;
        T_DoubleBuffer m_faceWeightsB;
        T_DoubleBuffer m_samplingCdf;
        T_DoubleBuffer m_samplingTotalWeight;
        T_DoubleBuffer m_systematicOffset;
        T_UnsignedBuffer m_stratifiedRayCounts;
        T_UnsignedBuffer m_stratifiedRayOffsets;
        T_UnsignedBuffer m_stratifiedRayFaces;
        T_CharBuffer m_samplingCdfScanBuffer;
        T_CharBuffer m_stratifiedCountScanBuffer;
        T_Reservoir m_reservoirA;
        T_Reservoir m_reservoirB;
        T_SamplingCdf m_samplingCdfSpans;
        bool m_inputReservoirA = true;
    };
} // namespace hase::core
