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
#include <limits>
#include <filesystem>
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
                   unsigned spectral)
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
            py::arg("spectral") = 0u)
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
        .def(
            "__repr__",
            [](ExperimentParameters const& p)
            {
                return "<ExperimentParameters minRaysPerSample=" + std::to_string(p.minRaysPerSample)
                       + ", maxRaysPerSample=" + std::to_string(p.maxRaysPerSample)
                       + ", mseThreshold=" + std::to_string(p.mseThreshold)
                       + ", useReflections=" + std::string(p.useReflections ? "True" : "False")
                       + ", spectral=" + std::to_string(p.spectral) + ">";
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
                   std::string deviceMode,
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
                    p.deviceMode = std::move(deviceMode);
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
            py::arg("deviceMode") = std::string("gpu"),
            py::arg("parallelMode") = std::string("threaded"),
            py::arg("writeVtk") = false,
            py::arg("devices") = std::vector<unsigned>{},
            py::arg("minSampleRange") = 0u,
            py::arg("maxSampleRange") = std::numeric_limits<unsigned>::max())
        .def_readwrite("maxRepetitions", &ComputeParameters::maxRepetitions)
        .def_readwrite("adaptiveSteps", &ComputeParameters::adaptiveSteps)
        .def_readwrite("maxGpus", &ComputeParameters::maxGpus)
        .def_readwrite("gpu_i", &ComputeParameters::gpu_i)
        .def_readwrite("deviceMode", &ComputeParameters::deviceMode)
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
                       + ", deviceMode='" + p.deviceMode + "', parallelMode='" + p.parallelMode + "'>";
            });

    py::class_<HostMesh>(m, "HostMesh")
        .def(
            py::init(
                [](std::vector<unsigned> triangleIndices,
                   unsigned numberOfTriangles,
                   unsigned numberOfLevels,
                   unsigned numberOfPoints,
                   float thicknessOfPrism,
                   std::vector<double> pointsVector,
                   std::vector<double> xOfTriangleCenter,
                   std::vector<double> yOfTriangleCenter,
                   std::vector<unsigned> positionsOfNormalVectors,
                   std::vector<double> xOfNormals,
                   std::vector<double> yOfNormals,
                   std::vector<int> forbiddenVector,
                   std::vector<int> neighborsVector,
                   std::vector<float> surfacesVector,
                   std::vector<double> betaValuesVector,
                   std::vector<double> betaCells,
                   std::vector<unsigned> cellTypes,
                   std::vector<float> refractiveIndices,
                   std::vector<float> reflectivities,
                   float nTot,
                   float crystalTFluo,
                   unsigned claddingNumber,
                   double claddingAbsorption)
                {
                    return HostMesh(
                        std::move(triangleIndices),
                        numberOfTriangles,
                        numberOfLevels,
                        numberOfPoints,
                        thicknessOfPrism,
                        std::move(pointsVector),
                        std::move(xOfTriangleCenter),
                        std::move(yOfTriangleCenter),
                        std::move(positionsOfNormalVectors),
                        std::move(xOfNormals),
                        std::move(yOfNormals),
                        std::move(forbiddenVector),
                        std::move(neighborsVector),
                        std::move(surfacesVector),
                        std::move(betaValuesVector),
                        std::move(betaCells),
                        std::move(cellTypes),
                        std::move(refractiveIndices),
                        std::move(reflectivities),
                        nTot,
                        crystalTFluo,
                        claddingNumber,
                        claddingAbsorption);
                }),
            py::arg("triangleIndices") = std::vector<unsigned>{},
            py::arg("numberOfTriangles") = 0u,
            py::arg("numberOfLevels") = 0u,
            py::arg("numberOfPoints") = 0u,
            py::arg("thicknessOfPrism") = 0.0f,
            py::arg("pointsVector") = std::vector<double>{},
            py::arg("xOfTriangleCenter") = std::vector<double>{},
            py::arg("yOfTriangleCenter") = std::vector<double>{},
            py::arg("positionsOfNormalVectors") = std::vector<unsigned>{},
            py::arg("xOfNormals") = std::vector<double>{},
            py::arg("yOfNormals") = std::vector<double>{},
            py::arg("forbiddenVector") = std::vector<int>{},
            py::arg("neighborsVector") = std::vector<int>{},
            py::arg("surfacesVector") = std::vector<float>{},
            py::arg("betaValuesVector") = std::vector<double>{},
            py::arg("betaCells") = std::vector<double>{},
            py::arg("cellTypes") = std::vector<unsigned>{},
            py::arg("refractiveIndices") = std::vector<float>{},
            py::arg("reflectivities") = std::vector<float>{},
            py::arg("nTot") = 0.0f,
            py::arg("crystalTFluo") = 0.0f,
            py::arg("claddingNumber") = 0u,
            py::arg("claddingAbsorption") = 0.0)
        .def_property_readonly("triangleIndices", [](HostMesh const& m) { return m.triangleIndices; })
        .def_property_readonly("numberOfTriangles", [](HostMesh const& m) { return m.numberOfTriangles; })
        .def_property_readonly("numberOfLevels", [](HostMesh const& m) { return m.numberOfLevels; })
        .def_property_readonly("numberOfPoints", [](HostMesh const& m) { return m.numberOfPoints; })
        .def_property_readonly("thicknessOfPrism", [](HostMesh const& m) { return m.thicknessOfPrism; })
        .def_readwrite("pointsVector", &HostMesh::pointsVector)
        .def_readwrite("xOfTriangleCenter", &HostMesh::xOfTriangleCenter)
        .def_readwrite("yOfTriangleCenter", &HostMesh::yOfTriangleCenter)
        .def_readwrite("positionsOfNormalVectors", &HostMesh::positionsOfNormalVectors)
        .def_readwrite("xOfNormals", &HostMesh::xOfNormals)
        .def_readwrite("yOfNormals", &HostMesh::yOfNormals)
        .def_readwrite("forbiddenVector", &HostMesh::forbiddenVector)
        .def_readwrite("neighborsVector", &HostMesh::neighborsVector)
        .def_readwrite("surfacesVector", &HostMesh::surfacesVector)
        .def_readwrite("betaValuesVector", &HostMesh::betaValuesVector)
        .def_readwrite("betaCells", &HostMesh::betaCells)
        .def_readwrite("cellTypes", &HostMesh::cellTypes)
        .def_readwrite("refractiveIndices", &HostMesh::refractiveIndices)
        .def_readwrite("reflectivities", &HostMesh::reflectivities)
        .def_readwrite("nTot", &HostMesh::nTot)
        .def_readwrite("crystalTFluo", &HostMesh::crystalTFluo)
        .def_readwrite("claddingNumber", &HostMesh::claddingNumber)
        .def_readwrite("claddingAbsorption", &HostMesh::claddingAbsorption)
        .def("toMesh", &HostMesh::toMesh)
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
        [](ExperimentParameters experiment, ComputeParameters compute, HostMesh host_mesh)
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
