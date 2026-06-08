/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <catch2/catch_get_random_seed.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <hase/hase.hpp>
#include <random/random.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <string>
#include <vector>
namespace fs = std::filesystem;
using TestApis
    = std::decay_t<decltype(alpaka::onHost::allBackends(alpaka::onHost::enabledApis, alpaka::exec::enabledExecutors))>;

#ifdef HASE_WORKDIR
static std::string workingDIR = HASE_WORKDIR;
#else
static std::string workingDIR = "..";
#endif
struct TestData
{
    std::string fileName;
    hase::core::ExperimentParameters expParams;
    hase::core::ComputeParameters compParams;
    std::optional<hase::core::HostMesh> mesh;
    hase::core::Result results;

    void createFromFile(fs::path const& cfg, fs::path const& dataFolder)
    {
        auto optionsMap = hase::parse::constructOptionsMapFromFile(cfg);

        optionsMap.set(hase::core::ExpSwitch::input_path, dataFolder.string());

        hase::parse::constructSimulationContextFromOptionsMap(optionsMap, expParams, compParams, mesh, results);
    }
};

auto makeTestData(fs::path const& cfg, fs::path const& dataFolder)
{
    std::vector<fs::path> inputFolders;

    for(auto const& entry : fs::directory_iterator(dataFolder))
    {
        if(entry.is_directory())
        {
            inputFolders.push_back(entry.path());
        }
    }

    std::vector<TestData> testDataEntries;
    testDataEntries.reserve(inputFolders.size());

    for(auto const& inputFolder : inputFolders)
    {
        TestData data;
        data.createFromFile(cfg, inputFolder);
        data.fileName = inputFolder.string();
        testDataEntries.push_back(std::move(data));
    }

    return testDataEntries;
}

void runSerial(TestData& testData)
{
    auto startHandle = hase::core::BaseVersionSerial{testData.expParams, testData.mesh.value(), testData.results};
    startHandle(testData.compParams.minSampleRange, testData.compParams.maxSampleRange);
}

void resetResultForSimulation(hase::core::Result& result)
{
    std::ranges::fill(result.phiAse, 0.0f);
    std::ranges::fill(result.mse, std::numeric_limits<double>::max());
    std::ranges::fill(result.totalRays, 0u);
    std::ranges::fill(result.dndtAse, 0.0);
}

fs::path comparisonCsvPath()
{
    return fs::path{workingDIR} / "tests/data/compSerial_Itest_results.csv";
}

void initializeComparisonCsv()
{
    std::ofstream csv{comparisonCsvPath()};
    csv << "dataset,backend,vector,index,serial,parallel,absDiff,relDiff\n";
}

std::vector<std::string> parseCsvLine(std::string const& line)
{
    std::vector<std::string> fields;
    std::string field;
    bool inQuotes = false;

    for(std::size_t i = 0; i < line.size(); ++i)
    {
        char const ch = line[i];
        if(inQuotes)
        {
            if(ch == '"')
            {
                if(i + 1 < line.size() && line[i + 1] == '"')
                {
                    field.push_back('"');
                    ++i;
                }
                else
                {
                    inQuotes = false;
                }
            }
            else
            {
                field.push_back(ch);
            }
        }
        else if(ch == '"')
        {
            inQuotes = true;
        }
        else if(ch == ',')
        {
            fields.push_back(field);
            field.clear();
        }
        else
        {
            field.push_back(ch);
        }
    }

    fields.push_back(field);
    return fields;
}

bool fillCachedSerialVector(
    std::string const& dataset,
    std::string const& vectorName,
    std::vector<double>& values,
    std::vector<bool>& seen,
    std::vector<std::string> const& fields)
{
    if(fields.size() != 8 || fields[0] != dataset || fields[2] != vectorName)
        return false;

    auto const index = static_cast<std::size_t>(std::stoull(fields[3]));
    if(index >= values.size())
        return false;

    values[index] = std::stod(fields[4]);
    seen[index] = true;
    return true;
}

bool fillCachedSerialVector(
    std::string const& dataset,
    std::string const& vectorName,
    std::vector<float>& values,
    std::vector<bool>& seen,
    std::vector<std::string> const& fields)
{
    if(fields.size() != 8 || fields[0] != dataset || fields[2] != vectorName)
        return false;

    auto const index = static_cast<std::size_t>(std::stoull(fields[3]));
    if(index >= values.size())
        return false;

    values[index] = std::stof(fields[4]);
    seen[index] = true;
    return true;
}

bool allCachedValuesPresent(std::vector<bool> const& seen, std::pair<unsigned, unsigned> const& sampleRange)
{
    for(std::size_t i = sampleRange.first; i <= sampleRange.second && i < seen.size(); ++i)
    {
        if(!seen[i])
            return false;
    }
    return true;
}

bool loadSerialResultsFromCsv(TestData& testData)
{
    auto const csvPath = comparisonCsvPath();
    if(!fs::exists(csvPath))
        return false;

    std::ifstream csv{csvPath};
    if(!csv)
        return false;

    auto cachedResults = testData.results;
    std::vector<bool> seenDndtAse(cachedResults.dndtAse.size(), false);
    std::vector<bool> seenPhiAse(cachedResults.phiAse.size(), false);

    std::string line;
    std::getline(csv, line);
    while(std::getline(csv, line))
    {
        if(line.empty())
            continue;

        auto const fields = parseCsvLine(line);
        fillCachedSerialVector(testData.fileName, "dndtAse", cachedResults.dndtAse, seenDndtAse, fields);
        fillCachedSerialVector(testData.fileName, "phiAse", cachedResults.phiAse, seenPhiAse, fields);
    }

    std::pair<unsigned, unsigned> const sampleRange{
        testData.compParams.minSampleRange,
        testData.compParams.maxSampleRange};
    if(!allCachedValuesPresent(seenDndtAse, sampleRange) || !allCachedValuesPresent(seenPhiAse, sampleRange))
        return false;

    testData.results = std::move(cachedResults);
    return true;
}

std::string csvQuote(std::string value)
{
    std::string quoted;
    quoted.reserve(value.size() + 2);
    quoted.push_back('"');
    for(char ch : value)
    {
        if(ch == '"')
        {
            quoted.push_back('"');
        }
        quoted.push_back(ch);
    }
    quoted.push_back('"');
    return quoted;
}

template<typename T>
void appendVectorComparisonCsv(
    std::string const& dataset,
    std::string const& backend,
    std::string const& vectorName,
    std::vector<T> const& expected,
    std::vector<T> const& actual,
    std::pair<unsigned, unsigned> const& sampleRange)
{
    if(expected.size() != actual.size())
        throw std::runtime_error("Vector size mismatch.");

    std::ofstream csv{comparisonCsvPath(), std::ios::app};
    csv << std::setprecision(17);

    for(std::size_t i = sampleRange.first; i <= sampleRange.second && i < expected.size(); ++i)
    {
        auto const e = static_cast<double>(expected[i]);
        auto const a = static_cast<double>(actual[i]);
        double const absDiff = std::abs(a - e);
        double const relDiff = absDiff / std::max(std::abs(e), 1e-30);

        csv << csvQuote(dataset) << ',' << csvQuote(backend) << ',' << csvQuote(vectorName) << ',' << i << ',' << e
            << ',' << a << ',' << absDiff << ',' << relDiff << '\n';
    }
}

template<typename T>
void assertVectorEqualWithin(
    std::string const& backend,
    std::string const& dataset,
    std::vector<T> const& expected,
    std::vector<T> const& actual,
    std::pair<unsigned, unsigned> const& sampleRange,
    double precision = 1e-2)
{
    if(expected.size() != actual.size())
        throw std::runtime_error("Vector size mismatch.");

    double sumExp = 0.0;
    double sumAct = 0.0;
    double sqDiff = 0.0;
    double sqExp = 0.0;
    double maxDiff = 0.0;

    for(std::size_t i = sampleRange.first; i <= sampleRange.second && i < expected.size(); ++i)
    {
        auto const e = static_cast<double>(expected[i]);
        auto const a = static_cast<double>(actual[i]);
        double const d = a - e;

        INFO("Non-finite vector comparison value");
        CAPTURE(backend);
        CAPTURE(dataset);
        CAPTURE(i);
        CAPTURE(e);
        CAPTURE(a);
        REQUIRE(std::isfinite(e));
        REQUIRE(std::isfinite(a));

        sumExp += e;
        sumAct += a;
        sqDiff += d * d;
        sqExp += e * e;
        maxDiff = std::max(maxDiff, std::abs(d));
    }

    double const relL2 = std::sqrt(sqDiff) / std::max(std::sqrt(sqExp), 1e-30);
    double const relSum = std::abs(sumAct - sumExp) / std::max(std::abs(sumExp), 1e-30);

    INFO("Relative Sum d comparison failed");
    CAPTURE(backend);
    CAPTURE(dataset);
    CAPTURE(maxDiff);
    CAPTURE(relL2);
    CAPTURE(relSum);
    CHECK(relSum < precision);
}

void compareResults(
    hase::core::Result const& expected,
    hase::core::Result const& actual,
    std::pair<unsigned, unsigned> sampleRange,
    std::string const& backend,
    std::string const& dataset)
{
    appendVectorComparisonCsv(dataset, backend, "dndtAse", expected.dndtAse, actual.dndtAse, sampleRange);
    appendVectorComparisonCsv(dataset, backend, "phiAse", expected.phiAse, actual.phiAse, sampleRange);

    assertVectorEqualWithin(backend, dataset, expected.dndtAse, actual.dndtAse, sampleRange);
    assertVectorEqualWithin(backend, dataset, expected.phiAse, actual.phiAse, sampleRange);
}

template<typename TBackend>
void runForEachBackend(hase::core::Result const& serialResults, TestData& currentData, TBackend&& backend)
{
    auto devSelector = alpaka::onHost::makeDeviceSelector(backend[alpaka::object::deviceSpec]);
    if(devSelector.getDeviceCount() == 0)
    {
        hase::core::dout(V_WARNING) << " No device found for backend." << std::endl;
        return;
    }
    auto sampleDevice = devSelector.makeDevice(0);
    currentData.compParams.backend = hase::core::getNameForBackend(backend, sampleDevice);
    std::cout << "running simulation with backend: " << currentData.compParams.backend << std::endl;
    resetResultForSimulation(currentData.results);
    hase::core::startSimulation<false>(
        currentData.expParams,
        currentData.compParams,
        currentData.results,
        currentData.mesh.value());
    std::cout << "simulation finished with backend: " << currentData.compParams.backend << std::endl;

    compareResults(
        serialResults,
        currentData.results,
        std::pair(currentData.compParams.minSampleRange, currentData.compParams.maxSampleRange),
        currentData.compParams.backend,
        currentData.fileName);
}

TEMPLATE_LIST_TEST_CASE("Test compare to Serial", "[compSerial]", TestApis)
{
    auto catch2Seed = Catch::getSeed();
    hase::random::SeedGenerator::get().updateSeed(catch2Seed);
    auto configPath = workingDIR + "/tests/data/cfg/compSerial.cfg";
    auto inputDataPath = workingDIR + "/example/c_example/input";
    auto testData = makeTestData(fs::path{configPath}, fs::path{inputDataPath});
    if(!fs::exists(comparisonCsvPath()))
    {
        initializeComparisonCsv();
    }

    auto backend = TestType::makeDict();

    for(auto& data : testData)
    {
        REQUIRE(data.expParams.useReflections == false);
        REQUIRE(data.expParams.maxRaysPerSample == data.expParams.minRaysPerSample);
        REQUIRE(data.expParams.monochromatic);
        REQUIRE(data.expParams.sigmaA.size() == 1);
        REQUIRE(data.expParams.sigmaE.size() == 1);
        if(loadSerialResultsFromCsv(data))
        {
            std::cout << "loaded serial results from: " << comparisonCsvPath() << std::endl;
        }
        else
        {
            runSerial(data);
        }
        //  Hard copy of the serial result before overwriting data.results
        hase::core::Result serialResults = data.results;
        runForEachBackend(serialResults, data, backend);
    }
}
