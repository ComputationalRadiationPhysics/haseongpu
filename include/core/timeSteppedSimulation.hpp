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
#include <core/simulation.hpp>
#include <core/simulationRunControl.hpp>
#include <core/simulationSnapshot.hpp>
#include <core/types.hpp>
#include <kernels/activePointMasks.hpp>
#include <kernels/derivativeComposition.hpp>
#include <kernels/oneDimensionalPump.hpp>
#include <kernels/pointBetaMapping.hpp>
#include <kernels/timeIntegrationUpdateKernels.hpp>

#include <algorithm>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace hase::core
{
    namespace detail
    {
        template<typename T_Buffer>
        std::vector<typename T_Buffer::value_type> copyToVector(auto const& queue, T_Buffer const& buffer)
        {
            auto hostBuffer = alpaka::onHost::allocHostLike(buffer);
            alpaka::onHost::memcpy(queue, hostBuffer, buffer);
            alpaka::onHost::wait(queue);
            auto const* data = alpaka::onHost::data(hostBuffer);
            return {data, data + buffer.getExtents().product()};
        }

        template<typename T_Buffer, typename T>
        void copyVectorToBuffer(auto const& queue, std::vector<T> const& values, T_Buffer& buffer)
        {
            auto hostView = alpaka::makeView(alpaka::api::host, values.data(), alpaka::Vec{values.size()});
            alpaka::onHost::memcpy(queue, buffer, hostView);
            alpaka::onHost::wait(queue);
        }

        double absorptionAtEmissionPeak(ExperimentParameters const& experiment)
        {
            if(experiment.sigmaA.empty() || experiment.sigmaE.empty())
            {
                return experiment.maxSigmaA;
            }
            auto const maxEmission = std::ranges::max_element(experiment.sigmaE);
            auto const index = static_cast<std::size_t>(std::distance(experiment.sigmaE.begin(), maxEmission));
            return experiment.sigmaA.at(std::min(index, experiment.sigmaA.size() - 1u));
        }

        void validateRunParameters(SimulationRunControl const& run)
        {
            if(run.timeStep <= 0.0)
            {
                throw std::runtime_error("simulation time_step must be positive");
            }
            if(run.numberOfSteps == 0u)
            {
                throw std::runtime_error("simulation number_of_steps must be positive");
            }
            if(run.pump.routine != PumpRoutine::ONE_DIMENSIONAL_Z_TRAVERSAL)
            {
                throw std::runtime_error("unsupported pump routine '" + run.pump.routine + "'");
            }
            if(run.pump.substeps < 2u)
            {
                throw std::runtime_error("pump_substeps must be at least 2");
            }
            if(run.pump.duration <= 0.0)
            {
                throw std::runtime_error("pump_duration must be positive");
            }
            if(run.pump.radiusX <= 0.0 || run.pump.radiusY <= 0.0)
            {
                throw std::runtime_error("pump radii must be positive");
            }
        }
    } // namespace detail

    template<typename T_Device, typename T_Executor>
    class CompiledSimulationRunner
    {
        using T_Queue = ALPAKA_TYPEOF(std::declval<T_Device>().makeQueue(alpaka::queueKind::blocking));
        using T_DoubleBuffer = ALPAKA_TYPEOF(alpaka::onHost::alloc<double>(std::declval<T_Device&>(), std::size_t{1}));
        using T_FloatBuffer = ALPAKA_TYPEOF(alpaka::onHost::alloc<float>(std::declval<T_Device&>(), std::size_t{1}));
        using T_UnsignedBuffer
            = ALPAKA_TYPEOF(alpaka::onHost::alloc<unsigned>(std::declval<T_Device&>(), std::size_t{1}));

    public:
        CompiledSimulationRunner(
            T_Device& device,
            T_Executor const& executor,
            ExperimentParameters& experiment,
            ComputeParameters& compute,
            SimulationRunControl const& run,
            HostMesh& hostMesh)
            : m_device(device)
            , m_queue(device.makeQueue(alpaka::queueKind::blocking))
            , m_devBundle(device, executor)
            , m_experiment(experiment)
            , m_compute(compute)
            , m_run(run)
            , m_hostMesh(hostMesh)
            , m_meshContainer(hostMesh.toDevice(device))
            , m_mesh(m_meshContainer.toView())
            , m_beta(hase::alpakaUtils::toDevice(m_queue, hostMesh.betaCells))
            , m_betaNext(alpaka::onHost::alloc<double>(device, static_cast<std::size_t>(m_mesh.numberOfSamples)))
            , m_stage(alpaka::onHost::alloc<double>(device, static_cast<std::size_t>(m_mesh.numberOfSamples)))
            , m_pumpedBeta(alpaka::onHost::alloc<double>(device, static_cast<std::size_t>(m_mesh.numberOfSamples)))
            , m_betaVolume(alpaka::onHost::alloc<double>(device, static_cast<std::size_t>(m_mesh.numberOfPrisms)))
            , m_phiAse(alpaka::onHost::alloc<float>(device, static_cast<std::size_t>(m_mesh.numberOfSamples)))
            , m_activeMask(alpaka::onHost::alloc<unsigned>(device, static_cast<std::size_t>(m_mesh.numberOfPoints)))
            , m_derivative(alpaka::onHost::alloc<double>(device, static_cast<std::size_t>(m_mesh.numberOfSamples)))
            , m_dndtPump(alpaka::onHost::alloc<double>(device, static_cast<std::size_t>(m_mesh.numberOfSamples)))
            , m_dndtAse(alpaka::onHost::alloc<double>(device, static_cast<std::size_t>(m_mesh.numberOfSamples)))
            , m_k1(alpaka::onHost::alloc<double>(device, static_cast<std::size_t>(m_mesh.numberOfSamples)))
            , m_k2(alpaka::onHost::alloc<double>(device, static_cast<std::size_t>(m_mesh.numberOfSamples)))
            , m_k3(alpaka::onHost::alloc<double>(device, static_cast<std::size_t>(m_mesh.numberOfSamples)))
            , m_k4(alpaka::onHost::alloc<double>(device, static_cast<std::size_t>(m_mesh.numberOfSamples)))
            , m_pumpForward(
                  alpaka::onHost::alloc<double>(
                      device,
                      static_cast<std::size_t>(m_mesh.numberOfPoints * m_mesh.numberOfLevels)))
            , m_pumpBackward(
                  alpaka::onHost::alloc<double>(
                      device,
                      static_cast<std::size_t>(m_mesh.numberOfPoints * m_mesh.numberOfLevels)))
        {
            hase::kernels::enqueueBuildActivePointMask(m_devBundle, m_queue, m_mesh, m_activeMask);
            alpaka::onHost::wait(m_queue);
        }

        void run(std::function<void(SimulationSnapshot const&)> const& callback)
        {
            for(unsigned step = 0u; step < m_run.numberOfSteps; ++step)
            {
                bool const pumpEnabled = step < m_run.pumpSteps;
                advanceOneStep(pumpEnabled);
                callback(makeSnapshot(step + 1u));
                alpaka::onHost::memcpy(m_queue, m_beta, m_betaNext);
                alpaka::onHost::wait(m_queue);
            }
        }

    private:
        void evaluateDerivative(auto& beta, bool pumpEnabled, bool refreshAse = true)
        {
            hase::kernels::enqueueMapPointBetaToPrismBeta(m_devBundle, m_queue, m_mesh, beta, m_betaVolume);
            alpaka::onHost::wait(m_queue);

            m_hostMesh.betaCells = detail::copyToVector(m_queue, beta);
            m_hostMesh.betaVolume = detail::copyToVector(m_queue, m_betaVolume);
            if(refreshAse)
            {
                initializeResult(m_run.enableAse ? 100000.0 : 0.0);

                if(m_run.enableAse)
                {
                    int const result
                        = hase::core::startSimulation<false>(m_experiment, m_compute, m_lastAseResult, m_hostMesh);
                    if(result != 0)
                    {
                        throw std::runtime_error("ASE evaluation failed with return code " + std::to_string(result));
                    }
                }
                detail::copyVectorToBuffer(m_queue, m_lastAseResult.phiAse, m_phiAse);
            }

            if(pumpEnabled)
            {
                hase::kernels::enqueueOneDimensionalPump(
                    m_devBundle,
                    m_queue,
                    m_mesh,
                    m_run.pump,
                    beta,
                    m_pumpedBeta,
                    m_pumpForward,
                    m_pumpBackward);
            }

            hase::kernels::enqueueComposeDerivative(
                m_devBundle,
                m_queue,
                m_mesh,
                detail::absorptionAtEmissionPeak(m_experiment),
                m_experiment.maxSigmaE,
                std::max(static_cast<double>(m_hostMesh.crystalTFluo), std::numeric_limits<double>::min()),
                m_run.pump.duration,
                pumpEnabled,
                beta,
                m_pumpedBeta,
                m_phiAse,
                m_activeMask,
                m_dndtPump,
                m_dndtAse,
                m_derivative);
            alpaka::onHost::wait(m_queue);
        }

        void advanceOneStep(bool pumpEnabled)
        {
            auto const& method = m_run.timeIntegration.method;
            if(method == TimeIntegrator::EXPLICIT_EULER)
            {
                evaluateDerivative(m_beta, pumpEnabled);
                enqueueAddScaled(m_beta, m_derivative, m_betaNext, m_run.timeStep);
            }
            else if(method == TimeIntegrator::HEUN)
            {
                evaluateDerivative(m_beta, pumpEnabled);
                alpaka::onHost::memcpy(m_queue, m_k1, m_derivative);
                enqueueAddScaled(m_beta, m_k1, m_stage, m_run.timeStep);
                evaluateDerivative(m_stage, pumpEnabled);
                enqueueHeun(m_beta, m_k1, m_derivative, m_betaNext);
            }
            else if(method == TimeIntegrator::MIDPOINT)
            {
                evaluateDerivative(m_beta, pumpEnabled);
                enqueueAddScaled(m_beta, m_derivative, m_stage, 0.5 * m_run.timeStep);
                evaluateDerivative(m_stage, pumpEnabled);
                enqueueAddScaled(m_beta, m_derivative, m_betaNext, m_run.timeStep);
            }
            else if(method == TimeIntegrator::RUNGE_KUTTA_4)
            {
                stepRungeKutta4(pumpEnabled);
            }
            else if(method == TimeIntegrator::FROZEN_PHI_ASE_RUNGE_KUTTA_4)
            {
                stepFrozenPhiAseRungeKutta4(pumpEnabled);
            }
            else if(method == TimeIntegrator::IMPLICIT_EULER)
            {
                stepImplicitEuler(pumpEnabled);
            }
            else if(method == TimeIntegrator::EXPONENTIAL_EULER)
            {
                evaluateDerivative(m_beta, pumpEnabled);
                enqueueExponentialEuler();
            }
            else
            {
                throw std::runtime_error("unsupported time integrator '" + method + "'");
            }

            enqueueClip(m_betaNext);
            hase::kernels::enqueueMapPointBetaToPrismBeta(m_devBundle, m_queue, m_mesh, m_betaNext, m_betaVolume);
            alpaka::onHost::wait(m_queue);
            m_hostMesh.betaCells = detail::copyToVector(m_queue, m_betaNext);
            m_hostMesh.betaVolume = detail::copyToVector(m_queue, m_betaVolume);
        }

        void stepRungeKutta4(bool pumpEnabled)
        {
            evaluateDerivative(m_beta, pumpEnabled);
            alpaka::onHost::memcpy(m_queue, m_k1, m_derivative);
            enqueueAddScaled(m_beta, m_k1, m_stage, 0.5 * m_run.timeStep);

            evaluateDerivative(m_stage, pumpEnabled);
            alpaka::onHost::memcpy(m_queue, m_k2, m_derivative);
            enqueueAddScaled(m_beta, m_k2, m_stage, 0.5 * m_run.timeStep);

            evaluateDerivative(m_stage, pumpEnabled);
            alpaka::onHost::memcpy(m_queue, m_k3, m_derivative);
            enqueueAddScaled(m_beta, m_k3, m_stage, m_run.timeStep);

            evaluateDerivative(m_stage, pumpEnabled);
            alpaka::onHost::memcpy(m_queue, m_k4, m_derivative);
            enqueueRungeKutta4();
        }

        void stepFrozenPhiAseRungeKutta4(bool pumpEnabled)
        {
            evaluateDerivative(m_beta, pumpEnabled, true);
            alpaka::onHost::memcpy(m_queue, m_k1, m_derivative);
            enqueueAddScaled(m_beta, m_k1, m_stage, 0.5 * m_run.timeStep);

            evaluateDerivative(m_stage, pumpEnabled, false);
            alpaka::onHost::memcpy(m_queue, m_k2, m_derivative);
            enqueueAddScaled(m_beta, m_k2, m_stage, 0.5 * m_run.timeStep);

            evaluateDerivative(m_stage, pumpEnabled, false);
            alpaka::onHost::memcpy(m_queue, m_k3, m_derivative);
            enqueueAddScaled(m_beta, m_k3, m_stage, m_run.timeStep);

            evaluateDerivative(m_stage, pumpEnabled, false);
            alpaka::onHost::memcpy(m_queue, m_k4, m_derivative);
            enqueueRungeKutta4();
        }

        void stepImplicitEuler(bool pumpEnabled)
        {
            alpaka::onHost::memcpy(m_queue, m_stage, m_beta);
            for(unsigned iteration = 0u; iteration < std::max(1u, m_run.timeIntegration.implicitIterations);
                ++iteration)
            {
                evaluateDerivative(m_stage, pumpEnabled);
                enqueueAddScaled(m_beta, m_derivative, m_betaNext, m_run.timeStep);
                alpaka::onHost::memcpy(m_queue, m_stage, m_betaNext);
            }
        }

        void initializeResult(double mseValue)
        {
            auto const numberOfSamples = m_hostMesh.numberOfPoints * m_hostMesh.numberOfLevels;
            m_lastAseResult = Result(
                std::vector<float>(numberOfSamples, 0.0f),
                std::vector<double>(numberOfSamples, mseValue),
                std::vector<unsigned>(numberOfSamples, 0u),
                std::vector<double>(numberOfSamples, 0.0));
        }

        SimulationSnapshot makeSnapshot(unsigned step)
        {
            auto dndtPump = detail::copyToVector(m_queue, m_dndtPump);
            auto dndtAse = detail::copyToVector(m_queue, m_dndtAse);
            m_lastAseResult.dndtAse = dndtAse;
            return SimulationSnapshot{
                step,
                static_cast<double>(step) * m_run.timeStep,
                m_hostMesh,
                m_lastAseResult,
                std::move(dndtPump),
                std::move(dndtAse)};
        }

        void enqueueAddScaled(auto& base, auto& slope, auto& out, double scale)
        {
            auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                m_devBundle.device,
                m_devBundle.executor,
                alpaka::Vec{m_mesh.numberOfSamples});
            m_queue.enqueue(
                frameSpec,
                alpaka::KernelBundle{hase::kernels::AddScaled{scale}, m_mesh, base, slope, out});
            alpaka::onHost::wait(m_queue);
        }

        void enqueueHeun(auto& base, auto& first, auto& second, auto& out)
        {
            auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                m_devBundle.device,
                m_devBundle.executor,
                alpaka::Vec{m_mesh.numberOfSamples});
            m_queue.enqueue(
                frameSpec,
                alpaka::KernelBundle{hase::kernels::CombineHeun{m_run.timeStep}, m_mesh, base, first, second, out});
            alpaka::onHost::wait(m_queue);
        }

        void enqueueRungeKutta4()
        {
            auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                m_devBundle.device,
                m_devBundle.executor,
                alpaka::Vec{m_mesh.numberOfSamples});
            m_queue.enqueue(
                frameSpec,
                alpaka::KernelBundle{
                    hase::kernels::CombineRungeKutta4{m_run.timeStep},
                    m_mesh,
                    m_beta,
                    m_k1,
                    m_k2,
                    m_k3,
                    m_k4,
                    m_betaNext});
            alpaka::onHost::wait(m_queue);
        }

        void enqueueExponentialEuler()
        {
            auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                m_devBundle.device,
                m_devBundle.executor,
                alpaka::Vec{m_mesh.numberOfSamples});
            m_queue.enqueue(
                frameSpec,
                alpaka::KernelBundle{
                    hase::kernels::ExponentialEulerUpdate{
                        m_run.timeStep,
                        std::max(static_cast<double>(m_hostMesh.crystalTFluo), std::numeric_limits<double>::min())},
                    m_mesh,
                    m_beta,
                    m_dndtPump,
                    m_dndtAse,
                    m_betaNext});
            alpaka::onHost::wait(m_queue);
        }

        void enqueueClip(auto& beta)
        {
            auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                m_devBundle.device,
                m_devBundle.executor,
                alpaka::Vec{m_mesh.numberOfSamples});
            m_queue.enqueue(frameSpec, alpaka::KernelBundle{hase::kernels::ClipBeta{}, m_mesh, beta});
            alpaka::onHost::wait(m_queue);
        }

        T_Device& m_device;
        T_Queue m_queue;
        hase::alpakaUtils::DevBundle<T_Device, T_Executor> m_devBundle;
        ExperimentParameters& m_experiment;
        ComputeParameters& m_compute;
        SimulationRunControl const& m_run;
        HostMesh& m_hostMesh;
        DeviceMeshContainer<T_Device> m_meshContainer;
        DeviceMeshView m_mesh;

        T_DoubleBuffer m_beta;
        T_DoubleBuffer m_betaNext;
        T_DoubleBuffer m_stage;
        T_DoubleBuffer m_pumpedBeta;
        T_DoubleBuffer m_betaVolume;
        T_FloatBuffer m_phiAse;
        T_UnsignedBuffer m_activeMask;
        T_DoubleBuffer m_derivative;
        T_DoubleBuffer m_dndtPump;
        T_DoubleBuffer m_dndtAse;
        T_DoubleBuffer m_k1;
        T_DoubleBuffer m_k2;
        T_DoubleBuffer m_k3;
        T_DoubleBuffer m_k4;
        T_DoubleBuffer m_pumpForward;
        T_DoubleBuffer m_pumpBackward;
        Result m_lastAseResult;
    };

    inline int startTimeSteppedSimulation(
        ExperimentParameters& experiment,
        ComputeParameters& compute,
        SimulationRunControl const& run,
        HostMesh& hostMesh,
        std::function<void(SimulationSnapshot const&)> const& callback)
    {
        detail::validateRunParameters(run);
        auto backends = alpaka::onHost::allBackends(alpaka::onHost::enabledApis, alpaka::exec::enabledExecutors);
        bool oneDidRun = false;
        alpaka::onHost::executeForEachIfHasDevice(
            [&](auto const& backend)
            {
                auto deviceSpec = backend[alpaka::object::deviceSpec];
                auto exec = backend[alpaka::object::exec];
                auto devSelector = alpaka::onHost::makeDeviceSelector(deviceSpec);
                if(devSelector.getDeviceCount() == 0u)
                {
                    return 0;
                }
                auto device = devSelector.makeDevice(0);
                if(!isSelected(backend, device, compute))
                {
                    return 0;
                }
                oneDidRun = true;
                CompiledSimulationRunner runner{device, exec, experiment, compute, run, hostMesh};
                runner.run(callback);
                return 0;
            },
            backends);

        return oneDidRun ? 0 : 1;
    }
} // namespace hase::core
