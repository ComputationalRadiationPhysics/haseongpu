#include <mesh.hpp>
#include <parser.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <types.hpp>

namespace py = pybind11;

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>
#include <simulation.hpp>

#include <filesystem>
#include <limits>
#include <string>
#include <utility>
#include <vector>

PYBIND11_MODULE(HASEonGPU, m)
{
    m.doc() = "Python bindings for HASEonGPU";

    py::class_<ExperimentParameters>(m, "ExperimentParameters")
        .def(
            py::init(
                [](unsigned minRaysPerSample,
                   unsigned maxRaysPerSample,
                   std::vector<double> lambdaA,
                   std::vector<double> lambdaE,
                   std::vector<double> sigmaA,
                   std::vector<double> sigmaE,
                   double maxSigmaA,
                   double maxSigmaE,
                   double mseThreshold,
                   bool useReflections,
                   unsigned spectral,
                   bool monochromatic)
                {
                    ExperimentParameters p{};
                    p.minRaysPerSample = minRaysPerSample;
                    p.maxRaysPerSample = maxRaysPerSample;
                    p.lambdaA = std::move(lambdaA);
                    p.lambdaE = std::move(lambdaE);
                    p.sigmaA = std::move(sigmaA);
                    p.sigmaE = std::move(sigmaE);
                    p.maxSigmaA = maxSigmaA;
                    p.maxSigmaE = maxSigmaE;
                    p.mseThreshold = mseThreshold;
                    p.useReflections = useReflections;
                    p.spectral = spectral;
                    p.monochromatic = monochromatic;
                    return p;
                }),
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
        .def_readwrite("minRaysPerSample", &ExperimentParameters::minRaysPerSample)
        .def_readwrite("maxRaysPerSample", &ExperimentParameters::maxRaysPerSample)
        .def_readwrite("lambdaA", &ExperimentParameters::lambdaA)
        .def_readwrite("lambdaE", &ExperimentParameters::lambdaE)
        .def_readwrite("sigmaA", &ExperimentParameters::sigmaA)
        .def_readwrite("sigmaE", &ExperimentParameters::sigmaE)
        .def_readwrite("maxSigmaA", &ExperimentParameters::maxSigmaA)
        .def_readwrite("maxSigmaE", &ExperimentParameters::maxSigmaE)
        .def_readwrite("mseThreshold", &ExperimentParameters::mseThreshold)
        .def_readwrite("useReflections", &ExperimentParameters::useReflections)
        .def_readwrite("spectral", &ExperimentParameters::spectral)
        .def_readwrite("monochromatic", &ExperimentParameters::monochromatic)
        .def(
            "__repr__",
            [](ExperimentParameters const& p)
            {
                return "<ExperimentParameters minRaysPerSample=" + std::to_string(p.minRaysPerSample)
                       + ", maxRaysPerSample=" + std::to_string(p.maxRaysPerSample)
                       + ", mseThreshold=" + std::to_string(p.mseThreshold) + ", useReflections="
                       + std::string(p.useReflections ? "True" : "False") + ", spectral=" + std::to_string(p.spectral)
                       + ", monochromatic=" + std::string(p.monochromatic ? "True" : "False") + ">";
            });

    py::class_<Result>(m, "Result")
        .def(
            py::init([](std::vector<float> phiAse,
                        std::vector<double> mse,
                        std::vector<unsigned> totalRays,
                        std::vector<double> dndtAse)
                     { return Result(std::move(phiAse), std::move(mse), std::move(totalRays), std::move(dndtAse)); }),
            py::arg("phiAse") = std::vector<float>{},
            py::arg("mse") = std::vector<double>{},
            py::arg("totalRays") = std::vector<unsigned>{},
            py::arg("dndtAse") = std::vector<double>{})
        .def_readwrite("phiAse", &Result::phiAse)
        .def_readwrite("mse", &Result::mse)
        .def_readwrite("totalRays", &Result::totalRays)
        .def_readwrite("dndtAse", &Result::dndtAse)
        .def_property_readonly("num_phiAse", [](Result const& r) { return r.phiAse.size(); })
        .def_property_readonly("num_mse", [](Result const& r) { return r.mse.size(); })
        .def_property_readonly("num_totalRays", [](Result const& r) { return r.totalRays.size(); })
        .def_property_readonly("num_dndtAse", [](Result const& r) { return r.dndtAse.size(); })
        .def(
            "__repr__",
            [](Result const& r)
            {
                return "<Result phiAse=" + std::to_string(r.phiAse.size()) + ", mse=" + std::to_string(r.mse.size())
                       + ", totalRays=" + std::to_string(r.totalRays.size())
                       + ", dndtAse=" + std::to_string(r.dndtAse.size()) + ">";
            });

    py::class_<ComputeParameters>(m, "ComputeParameters")
        .def(
            py::init(
                [](unsigned maxRepetitions,
                   unsigned adaptiveSteps,
                   unsigned maxGpus,
                   unsigned gpu_i,
                   std::string backend,
                   std::string parallelMode,
                   bool writeVtk,
                   std::vector<unsigned> devices,
                   unsigned minSampleRange,
                   unsigned maxSampleRange)
                {
                    ComputeParameters p{};
                    p.maxRepetitions = maxRepetitions;
                    p.adaptiveSteps = adaptiveSteps;
                    p.maxGpus = maxGpus;
                    p.gpu_i = gpu_i;
                    p.backend = std::move(backend);
                    p.parallelMode = std::move(parallelMode);
                    p.writeVtk = writeVtk;
                    p.devices = std::move(devices);
                    p.minSampleRange = minSampleRange;
                    p.maxSampleRange = maxSampleRange;
                    return p;
                }),
            py::arg("maxRepetitions") = 4u,
            py::arg("adaptiveSteps") = 4u,
            py::arg("maxGpus") = 1u,
            py::arg("gpu_i") = 0u,
            py::arg("backend") = std::string("gpu"),
            py::arg("parallelMode") = std::string("single"),
            py::arg("writeVtk") = false,
            py::arg("devices") = std::vector<unsigned>{},
            py::arg("minSampleRange") = 0u,
            py::arg("maxSampleRange") = std::numeric_limits<unsigned>::max())
        .def_readwrite("maxRepetitions", &ComputeParameters::maxRepetitions)
        .def_readwrite("adaptiveSteps", &ComputeParameters::adaptiveSteps)
        .def_readwrite("maxGpus", &ComputeParameters::maxGpus)
        .def_readwrite("gpu_i", &ComputeParameters::gpu_i)
        .def_readwrite("backend", &ComputeParameters::backend)
        .def_readwrite("parallelMode", &ComputeParameters::parallelMode)
        .def_readwrite("writeVtk", &ComputeParameters::writeVtk)
        .def_readwrite("devices", &ComputeParameters::devices)
        .def_readwrite("minSampleRange", &ComputeParameters::minSampleRange)
        .def_readwrite("maxSampleRange", &ComputeParameters::maxSampleRange)
        .def(
            "__repr__",
            [](ComputeParameters const& p)
            {
                return "<ComputeParameters maxRepetitions=" + std::to_string(p.maxRepetitions)
                       + ", adaptiveSteps=" + std::to_string(p.adaptiveSteps)
                       + ", maxGpus=" + std::to_string(p.maxGpus) + ", gpu_i=" + std::to_string(p.gpu_i)
                       + ", backend='" + p.backend + "', parallelMode='" + p.parallelMode + "'>";
            });

    py::class_<HostMesh>(m, "HostMesh")
        .def(
            py::init(
                [](std::vector<unsigned> trianglePointIndices,
                   unsigned numberOfTriangles,
                   unsigned numberOfLevels,
                   unsigned numberOfPoints,
                   float thickness,
                   std::vector<double> points,
                   std::vector<double> triangleCenterX,
                   std::vector<double> triangleCenterY,
                   std::vector<unsigned> triangleNormalPoint,
                   std::vector<double> triangleNormalsX,
                   std::vector<double> triangleNormalsY,
                   std::vector<int> forbiddenEdge,
                   std::vector<int> triangleNeighbors,
                   std::vector<float> triangleSurfaces,
                   std::vector<double> betaVolume,
                   std::vector<double> betaCells,
                   std::vector<unsigned> claddingCellTypes,
                   std::vector<float> refractiveIndices,
                   std::vector<float> reflectivities,
                   float nTot,
                   float crystalTFluo,
                   unsigned claddingNumber,
                   double claddingAbsorption)
                {
                    return HostMesh(
                        std::move(trianglePointIndices),
                        numberOfTriangles,
                        numberOfLevels,
                        numberOfPoints,
                        thickness,
                        std::move(points),
                        std::move(triangleCenterX),
                        std::move(triangleCenterY),
                        std::move(triangleNormalPoint),
                        std::move(triangleNormalsX),
                        std::move(triangleNormalsY),
                        std::move(forbiddenEdge),
                        std::move(triangleNeighbors),
                        std::move(triangleSurfaces),
                        std::move(betaVolume),
                        std::move(betaCells),
                        std::move(claddingCellTypes),
                        std::move(refractiveIndices),
                        std::move(reflectivities),
                        nTot,
                        crystalTFluo,
                        claddingNumber,
                        claddingAbsorption);
                }),
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
        .def_readwrite("trianglePointIndices", &HostMesh::trianglePointIndices)
        .def_readwrite("numberOfTriangles", &HostMesh::numberOfTriangles)
        .def_readwrite("numberOfLevels", &HostMesh::numberOfLevels)
        .def_readwrite("numberOfPoints", &HostMesh::numberOfPoints)
        .def_readwrite("thickness", &HostMesh::thickness)
        .def_readwrite("points", &HostMesh::points)
        .def_readwrite("triangleCenterX", &HostMesh::triangleCenterX)
        .def_readwrite("triangleCenterY", &HostMesh::triangleCenterY)
        .def_readwrite("triangleNormalPoint", &HostMesh::triangleNormalPoint)
        .def_readwrite("triangleNormalsX", &HostMesh::triangleNormalsX)
        .def_readwrite("triangleNormalsY", &HostMesh::triangleNormalsY)
        .def_readwrite("forbiddenEdge", &HostMesh::forbiddenEdge)
        .def_readwrite("triangleNeighbors", &HostMesh::triangleNeighbors)
        .def_readwrite("triangleSurfaces", &HostMesh::triangleSurfaces)
        .def_readwrite("betaVolume", &HostMesh::betaVolume)
        .def_readwrite("betaCells", &HostMesh::betaCells)
        .def_readwrite("claddingCellTypes", &HostMesh::claddingCellTypes)
        .def_readwrite("refractiveIndices", &HostMesh::refractiveIndices)
        .def_readwrite("reflectivities", &HostMesh::reflectivities)
        .def_readwrite("nTot", &HostMesh::nTot)
        .def_readwrite("crystalTFluo", &HostMesh::crystalTFluo)
        .def_readwrite("claddingNumber", &HostMesh::claddingNumber)
        .def_readwrite("claddingAbsorption", &HostMesh::claddingAbsorption)
        .def("calcTotalReflectionAngles", &HostMesh::calcTotalReflectionAngles)
        .def(
            "__repr__",
            [](HostMesh const& hm)
            {
                return "<HostMesh numberOfTriangles=" + std::to_string(hm.numberOfTriangles)
                       + ", numberOfLevels=" + std::to_string(hm.numberOfLevels)
                       + ", numberOfPoints=" + std::to_string(hm.numberOfPoints) + ">";
            });

    m.def(
        "calcPhiASE",
        [](ExperimentParameters& experiment, ComputeParameters& compute, HostMesh& host_mesh)
        {
            Result result;
            int const rc = pythonEntry(experiment, compute, result, host_mesh);
            if(rc != 0)
            {
                throw std::runtime_error("pythonEntry failed with return code " + std::to_string(rc));
            }
            return result;
        },
        py::arg("experiment"),
        py::arg("compute"),
        py::arg("host_mesh"));
}
