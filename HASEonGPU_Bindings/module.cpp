/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <core/cancellation.hpp>
#include <core/mesh.hpp>
#include <core/types.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <random/random.hpp>

namespace py = pybind11;

#include <core/simulation.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include <chrono>
#include <filesystem>
#include <future>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

PYBIND11_MODULE(HASEonGPU, m)
{
    m.doc() = "Python bindings for HASEonGPU";

    py::class_<hase::core::ExperimentParameters>(m, "ExperimentParameters")
        .def(
            py::init<
                unsigned,
                unsigned,
                std::vector<double>,
                std::vector<double>,
                std::vector<double>,
                std::vector<double>,
                double,
                double,
                double,
                bool,
                unsigned,
                bool>(),
            py::arg("minRaysPerSample") = 100000u,
            py::arg("maxRaysPerSample") = 100000u,
            py::arg("lambdaA") = std::vector<double>{},
            py::arg("lambdaE") = std::vector<double>{},
            py::arg("sigmaA") = std::vector<double>{},
            py::arg("sigmaE") = std::vector<double>{},
            py::arg("maxSigmaA") = 0.0,
            py::arg("maxSigmaE") = 0.0,
            py::arg("mseThreshold") = 0.1,
            py::arg("useReflections") = true,
            py::arg("spectral") = 0u,
            py::arg("monochromatic") = false)
        .def_readwrite("minRaysPerSample", &hase::core::ExperimentParameters::minRaysPerSample)
        .def_readwrite("maxRaysPerSample", &hase::core::ExperimentParameters::maxRaysPerSample)
        .def_readwrite("lambdaA", &hase::core::ExperimentParameters::lambdaA)
        .def_readwrite("lambdaE", &hase::core::ExperimentParameters::lambdaE)
        .def_readwrite("sigmaA", &hase::core::ExperimentParameters::sigmaA)
        .def_readwrite("sigmaE", &hase::core::ExperimentParameters::sigmaE)
        .def_readwrite("maxSigmaA", &hase::core::ExperimentParameters::maxSigmaA)
        .def_readwrite("maxSigmaE", &hase::core::ExperimentParameters::maxSigmaE)
        .def_readwrite("mseThreshold", &hase::core::ExperimentParameters::mseThreshold)
        .def_readwrite("useReflections", &hase::core::ExperimentParameters::useReflections)
        .def_readwrite("spectral", &hase::core::ExperimentParameters::spectral)
        .def_readwrite("monochromatic", &hase::core::ExperimentParameters::monochromatic)
        .def(
            "__repr__",
            [](hase::core::ExperimentParameters const& p)
            {
                return "<ExperimentParameters minRaysPerSample=" + std::to_string(p.minRaysPerSample)
                       + ", maxRaysPerSample=" + std::to_string(p.maxRaysPerSample)
                       + ", mseThreshold=" + std::to_string(p.mseThreshold) + ", useReflections="
                       + std::string(p.useReflections ? "True" : "False") + ", spectral=" + std::to_string(p.spectral)
                       + ", monochromatic=" + std::string(p.monochromatic ? "True" : "False") + ">";
            });

    py::class_<hase::core::Result>(m, "Result")
        .def(
            py::init<
                std::vector<float>,
                std::vector<double>,
                std::vector<unsigned>,
                std::vector<double>,
                std::vector<unsigned>>(),
            py::arg("phiAse") = std::vector<float>{},
            py::arg("mse") = std::vector<double>{},
            py::arg("totalRays") = std::vector<unsigned>{},
            py::arg("dndtAse") = std::vector<double>{},
            py::arg("droppedRays") = std::vector<unsigned>{})
        .def_readwrite("phiAse", &hase::core::Result::phiAse)
        .def_readwrite("mse", &hase::core::Result::mse)
        .def_readwrite("totalRays", &hase::core::Result::totalRays)
        .def_readwrite("dndtAse", &hase::core::Result::dndtAse)
        .def_readwrite("droppedRays", &hase::core::Result::droppedRays)
        .def_property_readonly("num_phiAse", [](hase::core::Result const& r) { return r.phiAse.size(); })
        .def_property_readonly("num_mse", [](hase::core::Result const& r) { return r.mse.size(); })
        .def_property_readonly("num_totalRays", [](hase::core::Result const& r) { return r.totalRays.size(); })
        .def_property_readonly("num_dndtAse", [](hase::core::Result const& r) { return r.dndtAse.size(); })
        .def_property_readonly("num_droppedRays", [](hase::core::Result const& r) { return r.droppedRays.size(); })
        .def(
            "__repr__",
            [](hase::core::Result const& r)
            {
                return "<Result phiAse=" + std::to_string(r.phiAse.size()) + ", mse=" + std::to_string(r.mse.size())
                       + ", totalRays=" + std::to_string(r.totalRays.size())
                       + ", dndtAse=" + std::to_string(r.dndtAse.size()) + ">";
            });

    py::class_<hase::core::ComputeParameters>(m, "ComputeParameters")
        .def(
            py::init(
                [](unsigned maxRepetitions,
                   unsigned adaptiveSteps,
                   unsigned numDevices,
                   std::string backend,
                   std::string parallelMode,
                   bool writeVtk,
                   std::vector<unsigned> devices,
                   unsigned minSampleRange,
                   unsigned maxSampleRange,
                   unsigned rngSeed)
                {
                    return hase::core::ComputeParameters(
                        maxRepetitions,
                        adaptiveSteps,
                        numDevices,
                        0u,
                        std::move(backend),
                        std::move(parallelMode),
                        writeVtk,
                        std::move(devices),
                        minSampleRange,
                        maxSampleRange,
                        rngSeed);
                }),
            py::arg("maxRepetitions") = 4u,
            py::arg("adaptiveSteps") = 4u,
            py::arg("numDevices") = 1u,
            py::arg("backend") = std::string("gpu"),
            py::arg("parallelMode") = std::string("single"),
            py::arg("writeVtk") = false,
            py::arg("devices") = std::vector<unsigned>{},
            py::arg("minSampleRange") = std::numeric_limits<unsigned>::max(),
            py::arg("maxSampleRange") = std::numeric_limits<unsigned>::max(),
            py::arg("rngSeed") = hase::core::ComputeParameters::unspecifiedRngSeed)
        .def_readwrite("maxRepetitions", &hase::core::ComputeParameters::maxRepetitions)
        .def_readwrite("adaptiveSteps", &hase::core::ComputeParameters::adaptiveSteps)
        .def_readwrite("numDevices", &hase::core::ComputeParameters::numDevices)
        .def_readwrite("backend", &hase::core::ComputeParameters::backend)
        .def_readwrite("parallelMode", &hase::core::ComputeParameters::parallelMode)
        .def_readwrite("writeVtk", &hase::core::ComputeParameters::writeVtk)
        .def_readwrite("devices", &hase::core::ComputeParameters::devices)
        .def_readwrite("minSampleRange", &hase::core::ComputeParameters::minSampleRange)
        .def_readwrite("maxSampleRange", &hase::core::ComputeParameters::maxSampleRange)
        .def_readwrite("rngSeed", &hase::core::ComputeParameters::rngSeed)
        .def(
            "__repr__",
            [](hase::core::ComputeParameters const& p)
            {
                return "<ComputeParameters maxRepetitions=" + std::to_string(p.maxRepetitions) + ", adaptiveSteps="
                       + std::to_string(p.adaptiveSteps) + ", numDevices=" + std::to_string(p.numDevices)
                       + ", backend='" + p.backend + "', parallelMode='" + p.parallelMode + "'>";
            });

    py::class_<hase::core::HostMesh>(m, "HostMesh")
        .def(
            py::init<
                std::vector<unsigned>,
                unsigned,
                unsigned,
                unsigned,
                float,
                std::vector<double>,
                std::vector<double>,
                std::vector<double>,
                std::vector<unsigned>,
                std::vector<double>,
                std::vector<double>,
                std::vector<int>,
                std::vector<int>,
                std::vector<float>,
                std::vector<double>,
                std::vector<double>,
                std::vector<unsigned>,
                std::vector<float>,
                std::vector<float>,
                float,
                float,
                unsigned,
                double>(),
            py::arg("trianglePointIndices") = std::vector<unsigned>{},
            py::arg("numberOfTriangles") = 0u,
            py::arg("numberOfLevels") = 0u,
            py::arg("numberOfPoints") = 0u,
            py::arg("thickness") = 0.0f,
            py::arg("points") = std::vector<double>{},
            py::arg("triangleCenterX") = std::vector<double>{},
            py::arg("triangleCenterY") = std::vector<double>{},
            py::arg("triangleNormalPoint") = std::vector<unsigned>{},
            py::arg("triangleNormalsX") = std::vector<double>{},
            py::arg("triangleNormalsY") = std::vector<double>{},
            py::arg("forbiddenEdge") = std::vector<int>{},
            py::arg("triangleNeighbors") = std::vector<int>{},
            py::arg("triangleSurfaces") = std::vector<float>{},
            py::arg("betaVolume") = std::vector<double>{},
            py::arg("betaCells") = std::vector<double>{},
            py::arg("claddingCellTypes") = std::vector<unsigned>{},
            py::arg("refractiveIndices") = std::vector<float>{},
            py::arg("reflectivities") = std::vector<float>{},
            py::arg("nTot") = 0.0f,
            py::arg("crystalTFluo") = 0.0f,
            py::arg("claddingNumber") = 0u,
            py::arg("claddingAbsorption") = 0.0)
        .def_readwrite("trianglePointIndices", &hase::core::HostMesh::trianglePointIndices)
        .def_readwrite("numberOfTriangles", &hase::core::HostMesh::numberOfTriangles)
        .def_readwrite("numberOfLevels", &hase::core::HostMesh::numberOfLevels)
        .def_readwrite("numberOfPoints", &hase::core::HostMesh::numberOfPoints)
        .def_readwrite("thickness", &hase::core::HostMesh::thickness)
        .def_readwrite("points", &hase::core::HostMesh::points)
        .def_readwrite("triangleCenterX", &hase::core::HostMesh::triangleCenterX)
        .def_readwrite("triangleCenterY", &hase::core::HostMesh::triangleCenterY)
        .def_readwrite("triangleNormalPoint", &hase::core::HostMesh::triangleNormalPoint)
        .def_readwrite("triangleNormalsX", &hase::core::HostMesh::triangleNormalsX)
        .def_readwrite("triangleNormalsY", &hase::core::HostMesh::triangleNormalsY)
        .def_readwrite("forbiddenEdge", &hase::core::HostMesh::forbiddenEdge)
        .def_readwrite("triangleNeighbors", &hase::core::HostMesh::triangleNeighbors)
        .def_readwrite("triangleSurfaces", &hase::core::HostMesh::triangleSurfaces)
        .def_readwrite("betaVolume", &hase::core::HostMesh::betaVolume)
        .def_readwrite("betaCells", &hase::core::HostMesh::betaCells)
        .def_readwrite("claddingCellTypes", &hase::core::HostMesh::claddingCellTypes)
        .def_readwrite("refractiveIndices", &hase::core::HostMesh::refractiveIndices)
        .def_readwrite("reflectivities", &hase::core::HostMesh::reflectivities)
        .def_readwrite("nTot", &hase::core::HostMesh::nTot)
        .def_readwrite("crystalTFluo", &hase::core::HostMesh::crystalTFluo)
        .def_readwrite("claddingNumber", &hase::core::HostMesh::claddingNumber)
        .def_readwrite("claddingAbsorption", &hase::core::HostMesh::claddingAbsorption)
        .def("calcTotalReflectionAngles", &hase::core::HostMesh::calcTotalReflectionAngles)
        .def(
            "__repr__",
            [](hase::core::HostMesh const& hm)
            {
                return "<HostMesh numberOfTriangles=" + std::to_string(hm.numberOfTriangles)
                       + ", numberOfLevels=" + std::to_string(hm.numberOfLevels)
                       + ", numberOfPoints=" + std::to_string(hm.numberOfPoints) + ">";
            });

    m.def(
        "calcPhiASE",
        [](hase::core::ExperimentParameters& experiment,
           hase::core::ComputeParameters& compute,
           hase::core::HostMesh& host_mesh)
        {
            hase::core::clearCancellation();
            auto future = std::async(
                std::launch::async,
                [experiment, compute, host_mesh]() mutable
                {
                    auto const numberOfSamples = host_mesh.numberOfPoints * host_mesh.numberOfLevels;
                    hase::core::Result result{
                        std::vector<float>(numberOfSamples, 0.0f),
                        std::vector<double>(numberOfSamples, 100000.0),
                        std::vector<unsigned>(numberOfSamples, 0u),
                        std::vector<double>(numberOfSamples, 0.0)};
                    int const rc = hase::core::startSimulation<false>(experiment, compute, result, host_mesh);
                    if(rc != 0)
                    {
                        throw std::runtime_error(
                            "hase::core::startSimulation failed with return code " + std::to_string(rc));
                    }
                    return result;
                });

            while(future.wait_for(std::chrono::milliseconds(100)) != std::future_status::ready)
            {
                if(PyErr_CheckSignals() != 0)
                {
                    hase::core::requestCancellation();
                    try
                    {
                        future.get();
                    }
                    catch(...)
                    {
                    }
                    hase::core::clearCancellation();
                    throw py::error_already_set();
                }
            }

            try
            {
                auto result = future.get();
                hase::core::clearCancellation();
                return result;
            }
            catch(...)
            {
                hase::core::clearCancellation();
                throw;
            }
        },
        py::arg("experiment"),
        py::arg("compute"),
        py::arg("host_mesh"));

    m.def("setRngSeed", [](unsigned seed) { hase::random::SeedGenerator::get().updateSeed(seed); }, py::arg("seed"));
    m.def("getRngSeed", []() { return hase::random::SeedGenerator::get().getSeed(); });
}
