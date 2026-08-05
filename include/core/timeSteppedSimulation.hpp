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
#include <alpakaUtils/backendNames.hpp>
#include <alpakaUtils/memory.hpp>
#include <alpakaUtils/utils.hpp>
#include <core/forwardPhiAseEvaluator.hpp>
#include <core/logging.hpp>
#include <core/mesh.hpp>
#include <core/simulationRunControl.hpp>
#include <core/simulationSnapshot.hpp>
#include <core/types.hpp>
#include <kernels/derivativeComposition.hpp>
#include <kernels/generalPump.hpp>
#include <kernels/timeIntegrationUpdateKernels.hpp>

#include <algorithm>
#ifdef HASE_ENABLE_STEP_TIMING
#    include <chrono>
#    include <cstdlib>
#    include <fstream>
#endif
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace hase::core
{
    namespace detail
    {
        inline double absorptionAtEmissionPeak(ExperimentParameters const& experiment)
        {
            if(experiment.sigmaA.empty() || experiment.sigmaE.empty())
            {
                return experiment.maxSigmaA;
            }
            auto const maxEmission = std::ranges::max_element(experiment.sigmaE);
            auto const index = static_cast<std::size_t>(std::distance(experiment.sigmaE.begin(), maxEmission));
            return experiment.sigmaA.at(std::min(index, experiment.sigmaA.size() - 1u));
        }

        inline std::vector<double> makeLumpedPointVolumes(HostMesh const& mesh)
        {
            std::vector<double> volumes(mesh.numberOfMeshPoints, 0.0);
            for(unsigned cell = 0u; cell < mesh.numberOfCells; ++cell)
            {
                double const share
                    = static_cast<double>(mesh.cellVolumes[cell]) / static_cast<double>(mesh.numberOfCellVertices);
                for(unsigned vertex = 0u; vertex < mesh.numberOfCellVertices; ++vertex)
                {
                    unsigned const point = mesh.cellPointIndices[cell * mesh.numberOfCellVertices + vertex];
                    volumes[point] += share;
                }
            }
            return volumes;
        }

        inline void validateRunParameters(SimulationRunControl const& run)
        {
            if(run.timeStep <= 0.0)
            {
                throw std::runtime_error("simulation time_step must be positive");
            }
            if(run.numberOfSteps == 0u)
            {
                throw std::runtime_error("simulation number_of_steps must be positive");
            }
            if(run.executionMode != SimulationExecutionMode::AUTONOMOUS
               && run.executionMode != SimulationExecutionMode::SYNCHRONIZED_DEBUG)
            {
                throw std::runtime_error("simulation execution_mode must be autonomous or synchronized-debug");
            }
            unsigned previousOutputStep = 0u;
            for(unsigned outputStep : run.outputSteps)
            {
                if(outputStep == 0u || outputStep > run.numberOfSteps)
                {
                    throw std::runtime_error(
                        "simulation output_steps entries must be completed-step indices in [1, number_of_steps]");
                }
                if(outputStep <= previousOutputStep)
                {
                    throw std::runtime_error("simulation output_steps must be strictly increasing and unique");
                }
                previousOutputStep = outputStep;
            }
            auto const supportedOutputFields = SimulationOutputField::all();
            if(run.outputFields.empty())
            {
                throw std::runtime_error("simulation output_fields must contain at least one field");
            }
            for(std::size_t index = 0u; index < run.outputFields.size(); ++index)
            {
                auto const& outputField = run.outputFields[index];
                if(std::ranges::find(supportedOutputFields, outputField) == supportedOutputFields.end())
                {
                    throw std::runtime_error("unsupported simulation output field '" + outputField + "'");
                }
                if(std::ranges::find(run.outputFields.begin(), run.outputFields.begin() + index, outputField)
                   != run.outputFields.begin() + index)
                {
                    throw std::runtime_error("simulation output_fields must be unique");
                }
            }
            if(run.executionMode == SimulationExecutionMode::AUTONOMOUS && !run.controlFields.empty())
                throw std::runtime_error("simulation control_fields require synchronized-debug mode");
            if(run.executionMode == SimulationExecutionMode::SYNCHRONIZED_DEBUG && !run.outputSteps.empty())
                throw std::runtime_error(
                    "synchronized-debug emits every completed step; output_steps must be omitted");
            auto const supportedControlFields = SimulationControlField::all();
            for(std::size_t index = 0u; index < run.controlFields.size(); ++index)
            {
                auto const& controlField = run.controlFields[index];
                if(std::ranges::find(supportedControlFields, controlField) == supportedControlFields.end())
                    throw std::runtime_error("unsupported simulation control field '" + controlField + "'");
                if(std::ranges::find(run.controlFields.begin(), run.controlFields.begin() + index, controlField)
                   != run.controlFields.begin() + index)
                    throw std::runtime_error("simulation control_fields must be unique");
            }
            if(run.pump.schemaVersion != 1u || run.pump.rayCount == 0u || run.pump.sources.empty())
                throw std::runtime_error("invalid general pump configuration");
            for(auto const& source : run.pump.sources)
                if(source.totalPower <= 0.0 || source.surfaces.empty() || source.wavelengths.empty()
                   || source.wavelengths.size() != source.spectralWeights.size()
                   || source.wavelengths.size() != source.sigmaAbsorption.size()
                   || source.wavelengths.size() != source.sigmaEmission.size() || source.polarAngles.empty()
                   || source.polarAngles.size() != source.azimuthalAngles.size()
                   || source.polarAngles.size() != source.angularWeights.size())
                    throw std::runtime_error("invalid general pump source configuration");
        }
    } // namespace detail

    template<typename T_Device, typename T_Executor>
    class CompiledSimulationRunner
    {
        using T_Queue = ALPAKA_TYPEOF(std::declval<T_Device>().makeQueue(alpaka::queueKind::blocking));
        using T_DoubleBuffer = ALPAKA_TYPEOF(alpaka::onHost::alloc<double>(std::declval<T_Device&>(), std::size_t{1}));
        using T_FloatBuffer = ALPAKA_TYPEOF(alpaka::onHost::alloc<float>(std::declval<T_Device&>(), std::size_t{1}));

    public:
        CompiledSimulationRunner(
            std::vector<T_Device> devices,
            T_Executor const& executor,
            ExperimentParameters& experiment,
            ComputeParameters& compute,
            SimulationRunControl const& run,
            HostMesh& hostMesh)
            : m_forwardAseContext(std::move(devices), executor, experiment, hostMesh)
            , m_device(m_forwardAseContext.primaryDevice())
            , m_queue(m_device.makeQueue(alpaka::queueKind::blocking))
            , m_devBundle(m_device, executor)
            , m_experiment(experiment)
            , m_compute(compute)
            , m_run(run)
            , m_hostMesh(hostMesh)
            , m_mesh(m_forwardAseContext.primaryMesh().toView())
            , m_beta(hase::alpakaUtils::toDevice(m_queue, hostMesh.betaVolume))
            , m_betaNext(alpaka::onHost::alloc<double>(m_device, static_cast<std::size_t>(m_mesh.numberOfCells)))
            , m_stage(alpaka::onHost::alloc<double>(m_device, static_cast<std::size_t>(m_mesh.numberOfCells)))
            , m_phiAse(alpaka::onHost::alloc<float>(m_device, static_cast<std::size_t>(m_mesh.numberOfCells)))
            , m_derivative(alpaka::onHost::alloc<double>(m_device, static_cast<std::size_t>(m_mesh.numberOfCells)))
            , m_dndtPump(alpaka::onHost::alloc<double>(m_device, static_cast<std::size_t>(m_mesh.numberOfCells)))
            , m_dndtAse(alpaka::onHost::alloc<double>(m_device, static_cast<std::size_t>(m_mesh.numberOfCells)))
            , m_k1(alpaka::onHost::alloc<double>(m_device, static_cast<std::size_t>(m_mesh.numberOfCells)))
            , m_k2(alpaka::onHost::alloc<double>(m_device, static_cast<std::size_t>(m_mesh.numberOfCells)))
            , m_k3(alpaka::onHost::alloc<double>(m_device, static_cast<std::size_t>(m_mesh.numberOfCells)))
            , m_k4(alpaka::onHost::alloc<double>(m_device, static_cast<std::size_t>(m_mesh.numberOfCells)))
            , m_cellPumpIntegral(
                  alpaka::onHost::alloc<double>(m_device, static_cast<std::size_t>(m_mesh.numberOfCells)))
            , m_pointPumpIntegral(
                  alpaka::onHost::alloc<double>(m_device, static_cast<std::size_t>(m_mesh.numberOfMeshPoints)))
            , m_lumpedPointVolume(hase::alpakaUtils::toDevice(m_queue, detail::makeLumpedPointVolumes(hostMesh)))
        {
            if(hostMesh.betaVolume.size() != hostMesh.numberOfCells)
                throw std::runtime_error("simulation beta_volume must contain exactly one value per cell");
        }

        void run(
            std::function<void(SimulationSnapshot const&)> const& callback,
            std::function<std::vector<double>(unsigned)> const& receiveControl)
        {
#ifdef HASE_ENABLE_STEP_TIMING
            std::ofstream timingCsv;
            if(auto const* timingPath = std::getenv("HASE_STEP_TIMING_CSV"))
            {
                timingCsv.open(timingPath);
                timingCsv << "revision,backend,step,elapsed_seconds,pump_enabled,ase_enabled\n";
            }
            char const* revision = std::getenv("HASE_BENCHMARK_REVISION");
            char const* backend = std::getenv("HASE_BENCHMARK_BACKEND");
#endif
            for(unsigned step = 0u; step < m_run.numberOfSteps; ++step)
            {
                bool const pumpEnabled = step < m_run.pumpSteps;
                bool const aseEnabled = m_run.enableAse && !(m_run.prePump && step == 0u);
#ifdef HASE_ENABLE_STEP_TIMING
                auto const started = std::chrono::steady_clock::now();
#endif
                advanceOneStep(pumpEnabled, aseEnabled);
#ifdef HASE_ENABLE_STEP_TIMING
                alpaka::onHost::wait(m_queue);
                if(timingCsv)
                {
                    std::chrono::duration<double> const elapsed = std::chrono::steady_clock::now() - started;
                    timingCsv << (revision ? revision : "") << ',' << (backend ? backend : "") << ',' << (step + 1u)
                              << ',' << elapsed.count() << ',' << (pumpEnabled ? 1 : 0) << ',' << (aseEnabled ? 1 : 0)
                              << '\n';
                    timingCsv.flush();
                }
#endif
                std::swap(m_beta, m_betaNext);
                unsigned const completedStep = step + 1u;
                bool const synchronizedDebug = m_run.executionMode == SimulationExecutionMode::SYNCHRONIZED_DEBUG;
                if(synchronizedDebug || shouldOutput(completedStep))
                {
                    callback(makeSnapshot(completedStep));
                }
                if(synchronizedDebug && completedStep < m_run.numberOfSteps)
                {
                    if(!receiveControl)
                        throw std::runtime_error("synchronized-debug requires a control receiver");
                    auto betaVolume = receiveControl(completedStep);
                    if(includesControl(SimulationControlField::BETA_VOLUME))
                    {
                        if(betaVolume.size() != m_mesh.numberOfCells)
                            throw std::runtime_error("synchronized beta_volume control has the wrong cell count");
                        detail::copyVectorToBuffer(m_queue, betaVolume, m_beta);
                    }
                }
            }
        }

    private:
        struct DerivativeBuffers
        {
            T_DoubleBuffer& betaVolume;
            T_FloatBuffer& phiAse;
            T_DoubleBuffer& dndtPump;
            T_DoubleBuffer& dndtAse;
            T_DoubleBuffer& derivative;
        };

        void evaluateDerivative(auto& beta, bool pumpEnabled, bool aseEnabled, bool refreshAse = true)
        {
            if(refreshAse)
            {
                initializeResult(aseEnabled ? 100000.0 : 0.0, m_hostMesh.numberOfCells);

                if(aseEnabled)
                {
                    m_phiAseDeviceResident
                        = m_forwardAseContext.evaluate(m_experiment, m_compute, m_hostMesh, beta, m_lastAseResult)
                              .deviceResidentPhi;
                    if(!m_phiAseDeviceResident)
                        detail::copyVectorToBuffer(m_queue, m_lastAseResult.phiAse, m_phiAse);
                }
                else
                {
                    m_phiAseDeviceResident = false;
                    alpaka::onHost::fill(
                        m_queue,
                        m_phiAse,
                        0.0f,
                        alpaka::Vec{static_cast<std::size_t>(m_mesh.numberOfCells)});
                }
            }

            if(pumpEnabled)
            {
                hase::kernels::enqueueGeneralPump(
                    m_devBundle,
                    m_queue,
                    m_hostMesh,
                    m_mesh,
                    m_run.pump,
                    beta,
                    m_cellPumpIntegral,
                    m_pointPumpIntegral,
                    m_lumpedPointVolume,
                    m_dndtPump);
            }

            auto& activePhiAse = m_phiAseDeviceResident ? m_forwardAseContext.primaryVolumePhiAse() : m_phiAse;
            DerivativeBuffers derivativeBuffers{beta, activePhiAse, m_dndtPump, m_dndtAse, m_derivative};
            hase::kernels::enqueueComposeDerivative(
                m_devBundle,
                m_queue,
                m_mesh,
                detail::absorptionAtEmissionPeak(m_experiment),
                m_experiment.maxSigmaE,
                std::max(static_cast<double>(m_hostMesh.crystalTFluo), std::numeric_limits<double>::min()),
                pumpEnabled,
                derivativeBuffers);
        }

        void advanceOneStep(bool pumpEnabled, bool aseEnabled)
        {
            auto const& method = m_run.timeIntegration.method;
            if(method == TimeIntegrator::EXPLICIT_EULER)
            {
                evaluateDerivative(m_beta, pumpEnabled, aseEnabled);
                enqueueAddScaled(m_beta, m_derivative, m_betaNext, m_run.timeStep);
            }
            else if(method == TimeIntegrator::HEUN)
            {
                evaluateDerivative(m_beta, pumpEnabled, aseEnabled);
                alpaka::onHost::memcpy(m_queue, m_k1, m_derivative);
                enqueueAddScaled(m_beta, m_k1, m_stage, m_run.timeStep);
                evaluateDerivative(m_stage, pumpEnabled, aseEnabled);
                enqueueHeun(m_beta, m_k1, m_derivative, m_betaNext);
            }
            else if(method == TimeIntegrator::MIDPOINT)
            {
                evaluateDerivative(m_beta, pumpEnabled, aseEnabled);
                enqueueAddScaled(m_beta, m_derivative, m_stage, 0.5 * m_run.timeStep);
                evaluateDerivative(m_stage, pumpEnabled, aseEnabled);
                enqueueAddScaled(m_beta, m_derivative, m_betaNext, m_run.timeStep);
            }
            else if(method == TimeIntegrator::RUNGE_KUTTA_4)
            {
                stepRungeKutta4(pumpEnabled, aseEnabled);
            }
            else if(method == TimeIntegrator::FROZEN_PHI_ASE_RUNGE_KUTTA_4)
            {
                stepFrozenPhiAseRungeKutta4(pumpEnabled, aseEnabled);
            }
            else if(method == TimeIntegrator::IMPLICIT_EULER)
            {
                stepImplicitEuler(pumpEnabled, aseEnabled);
            }
            else if(method == TimeIntegrator::EXPONENTIAL_EULER)
            {
                evaluateDerivative(m_beta, pumpEnabled, aseEnabled);
                enqueueExponentialEuler();
            }
            else
            {
                throw std::runtime_error("unsupported time integrator '" + method + "'");
            }

            enqueueClip(m_betaNext);
        }

        void stepRungeKutta4(bool pumpEnabled, bool aseEnabled)
        {
            evaluateDerivative(m_beta, pumpEnabled, aseEnabled);
            alpaka::onHost::memcpy(m_queue, m_k1, m_derivative);
            enqueueAddScaled(m_beta, m_k1, m_stage, 0.5 * m_run.timeStep);

            evaluateDerivative(m_stage, pumpEnabled, aseEnabled);
            alpaka::onHost::memcpy(m_queue, m_k2, m_derivative);
            enqueueAddScaled(m_beta, m_k2, m_stage, 0.5 * m_run.timeStep);

            evaluateDerivative(m_stage, pumpEnabled, aseEnabled);
            alpaka::onHost::memcpy(m_queue, m_k3, m_derivative);
            enqueueAddScaled(m_beta, m_k3, m_stage, m_run.timeStep);

            evaluateDerivative(m_stage, pumpEnabled, aseEnabled);
            alpaka::onHost::memcpy(m_queue, m_k4, m_derivative);
            enqueueRungeKutta4();
        }

        void stepFrozenPhiAseRungeKutta4(bool pumpEnabled, bool aseEnabled)
        {
            evaluateDerivative(m_beta, pumpEnabled, aseEnabled, true);
            alpaka::onHost::memcpy(m_queue, m_k1, m_derivative);
            enqueueAddScaled(m_beta, m_k1, m_stage, 0.5 * m_run.timeStep);

            evaluateDerivative(m_stage, pumpEnabled, aseEnabled, false);
            alpaka::onHost::memcpy(m_queue, m_k2, m_derivative);
            enqueueAddScaled(m_beta, m_k2, m_stage, 0.5 * m_run.timeStep);

            evaluateDerivative(m_stage, pumpEnabled, aseEnabled, false);
            alpaka::onHost::memcpy(m_queue, m_k3, m_derivative);
            enqueueAddScaled(m_beta, m_k3, m_stage, m_run.timeStep);

            evaluateDerivative(m_stage, pumpEnabled, aseEnabled, false);
            alpaka::onHost::memcpy(m_queue, m_k4, m_derivative);
            enqueueRungeKutta4();
        }

        void stepImplicitEuler(bool pumpEnabled, bool aseEnabled)
        {
            alpaka::onHost::memcpy(m_queue, m_stage, m_beta);
            for(unsigned iteration = 0u; iteration < std::max(1u, m_run.timeIntegration.implicitIterations);
                ++iteration)
            {
                evaluateDerivative(m_stage, pumpEnabled, aseEnabled);
                enqueueAddScaled(m_beta, m_derivative, m_betaNext, m_run.timeStep);
                alpaka::onHost::memcpy(m_queue, m_stage, m_betaNext);
            }
        }

        void initializeResult(double standardErrorValue, unsigned resultSize)
        {
            m_lastAseResult = Result(
                std::vector<float>(resultSize, 0.0f),
                std::vector<double>(resultSize, standardErrorValue),
                std::vector<double>(resultSize, 0.0),
                std::vector<unsigned>(resultSize, 0u),
                std::vector<double>(resultSize, 0.0));
        }

        [[nodiscard]] bool includes(std::string const& field) const
        {
            return std::ranges::find(m_run.outputFields, field) != m_run.outputFields.end();
        }

        [[nodiscard]] bool includesControl(std::string const& field) const
        {
            return std::ranges::find(m_run.controlFields, field) != m_run.controlFields.end();
        }

        [[nodiscard]] bool shouldOutput(unsigned completedStep)
        {
            if(m_run.outputSteps.empty())
            {
                return true;
            }
            if(m_nextOutputStep >= m_run.outputSteps.size() || m_run.outputSteps[m_nextOutputStep] != completedStep)
            {
                return false;
            }
            ++m_nextOutputStep;
            return true;
        }

        static void copyStatus(Result const& source, Result& target)
        {
            target.srmStatus = source.srmStatus;
            target.srmPasses = source.srmPasses;
            target.srmRemainingFraction = source.srmRemainingFraction;
            target.srmMaxIterations = source.srmMaxIterations;
            target.srmDivergenceStreak = source.srmDivergenceStreak;
        }

        SimulationSnapshot makeSnapshot(unsigned step)
        {
            std::vector<double> betaVolume;
            std::vector<double> dndtPump;
            std::vector<double> dndtAse;
            Result aseResult;
            copyStatus(m_lastAseResult, aseResult);

            if(includes(SimulationOutputField::BETA_VOLUME))
                betaVolume = detail::copyToVector(m_queue, m_beta);
            bool const includePhiAse = includes(SimulationOutputField::PHI_ASE);
            bool const includeStandardError = includes(SimulationOutputField::STANDARD_ERROR);
            bool const includeRelativeStandardError = includes(SimulationOutputField::RELATIVE_STANDARD_ERROR);
            bool const includeTotalRays = includes(SimulationOutputField::TOTAL_RAYS);
            if(m_phiAseDeviceResident)
            {
                auto deviceResult = m_forwardAseContext.downloadPrimaryResult(
                    includePhiAse,
                    includeStandardError,
                    includeRelativeStandardError,
                    includeTotalRays);
                aseResult.phiAse = std::move(deviceResult.phiAse);
                aseResult.standardError = std::move(deviceResult.standardError);
                aseResult.relativeStandardError = std::move(deviceResult.relativeStandardError);
                aseResult.totalRays = std::move(deviceResult.totalRays);
                aseResult.droppedRays = std::move(deviceResult.droppedRays);
            }
            else
            {
                if(includePhiAse)
                    aseResult.phiAse = m_lastAseResult.phiAse;
                if(includeStandardError)
                    aseResult.standardError = m_lastAseResult.standardError;
                if(includeRelativeStandardError)
                    aseResult.relativeStandardError = m_lastAseResult.relativeStandardError;
                if(includeTotalRays)
                {
                    aseResult.totalRays = m_lastAseResult.totalRays;
                    aseResult.droppedRays = m_lastAseResult.droppedRays;
                }
            }
            if(includes(SimulationOutputField::DNDT_ASE))
            {
                dndtAse = detail::copyToVector(m_queue, m_dndtAse);
                aseResult.dndtAse = dndtAse;
            }
            if(includes(SimulationOutputField::DNDT_PUMP))
                dndtPump = detail::copyToVector(m_queue, m_dndtPump);

            return SimulationSnapshot{
                step,
                static_cast<double>(step) * m_run.timeStep,
                std::move(betaVolume),
                std::move(aseResult),
                std::move(dndtPump),
                std::move(dndtAse),
                m_run.outputFields};
        }

        void enqueueAddScaled(auto& base, auto& slope, auto& out, double scale)
        {
            auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                m_devBundle.device,
                m_devBundle.executor,
                alpaka::Vec{m_mesh.numberOfCells});
            m_queue.enqueue(
                frameSpec,
                alpaka::KernelBundle{hase::kernels::AddScaled{scale}, m_mesh, base, slope, out});
        }

        void enqueueHeun(auto& base, auto& first, auto& second, auto& out)
        {
            auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                m_devBundle.device,
                m_devBundle.executor,
                alpaka::Vec{m_mesh.numberOfCells});
            m_queue.enqueue(
                frameSpec,
                alpaka::KernelBundle{hase::kernels::CombineHeun{m_run.timeStep}, m_mesh, base, first, second, out});
        }

        void enqueueRungeKutta4()
        {
            auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                m_devBundle.device,
                m_devBundle.executor,
                alpaka::Vec{m_mesh.numberOfCells});
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
        }

        void enqueueExponentialEuler()
        {
            auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                m_devBundle.device,
                m_devBundle.executor,
                alpaka::Vec{m_mesh.numberOfCells});
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
        }

        void enqueueClip(auto& beta)
        {
            auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                m_devBundle.device,
                m_devBundle.executor,
                alpaka::Vec{m_mesh.numberOfCells});
            m_queue.enqueue(frameSpec, alpaka::KernelBundle{hase::kernels::ClipBeta{}, m_mesh, beta});
        }

        ForwardPhiAseContext<T_Device, T_Executor> m_forwardAseContext;
        T_Device& m_device;
        T_Queue m_queue;
        hase::alpakaUtils::DevBundle<T_Device, T_Executor> m_devBundle;
        ExperimentParameters& m_experiment;
        ComputeParameters& m_compute;
        SimulationRunControl const& m_run;
        HostMesh& m_hostMesh;
        DeviceMeshView m_mesh;

        T_DoubleBuffer m_beta;
        T_DoubleBuffer m_betaNext;
        T_DoubleBuffer m_stage;
        T_FloatBuffer m_phiAse;
        T_DoubleBuffer m_derivative;
        T_DoubleBuffer m_dndtPump;
        T_DoubleBuffer m_dndtAse;
        T_DoubleBuffer m_k1;
        T_DoubleBuffer m_k2;
        T_DoubleBuffer m_k3;
        T_DoubleBuffer m_k4;
        T_DoubleBuffer m_cellPumpIntegral;
        T_DoubleBuffer m_pointPumpIntegral;
        T_DoubleBuffer m_lumpedPointVolume;
        Result m_lastAseResult;
        std::size_t m_nextOutputStep = 0u;
        bool m_phiAseDeviceResident = false;
    };

    inline int startTimeSteppedSimulation(
        ExperimentParameters& experiment,
        ComputeParameters& compute,
        SimulationRunControl const& run,
        HostMesh& hostMesh,
        std::function<void(SimulationSnapshot const&)> const& callback,
        std::function<std::vector<double>(unsigned)> const& receiveControl = {})
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
                auto sampleDevice = devSelector.makeDevice(0);
                if(hase::alpakaUtils::getNameForBackend(backend, sampleDevice) != compute.backend)
                {
                    return 0;
                }
                std::size_t const deviceCount = devSelector.getDeviceCount();
                compute.devices.resize(deviceCount);
                std::iota(compute.devices.begin(), compute.devices.end(), 0u);
                compute.gpu_i = compute.devices.front();
                if(compute.numDevices == 0u)
                    compute.numDevices = static_cast<unsigned>(deviceCount);
                if(compute.numDevices > deviceCount)
                {
                    dout(V_WARNING) << "Requested number of devices (" << compute.numDevices
                                    << ") exceeds the available device count (" << deviceCount
                                    << "); using all available devices." << std::endl;
                    compute.numDevices = static_cast<unsigned>(deviceCount);
                }
                compute.devices.resize(compute.numDevices);
                using T_Device = ALPAKA_TYPEOF(sampleDevice);
                std::vector<T_Device> devices;
                devices.reserve(compute.devices.size());
                for(unsigned deviceIndex : compute.devices)
                    devices.emplace_back(devSelector.makeDevice(deviceIndex));
                oneDidRun = true;
                CompiledSimulationRunner runner{std::move(devices), exec, experiment, compute, run, hostMesh};
                runner.run(callback, receiveControl);
                return 0;
            },
            backends);

        if(!oneDidRun)
        {
            std::ostringstream message;
            message << "Backend '" << compute.backend
                    << "' did not match any available backend with an available device. Available backends:";
            for(auto const& element : hase::alpakaUtils::availableBackendNames())
            {
                message << "\n  " << element;
            }
            throw std::runtime_error(message.str());
        }

        return 0;
    }
} // namespace hase::core
