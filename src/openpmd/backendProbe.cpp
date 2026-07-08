#include <openPMD/openPMD.hpp>
#include <openpmd/backendProbe.hpp>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

namespace io = openPMD;

namespace
{
    struct Backend
    {
        std::string_view name;
        std::string_view suffix;
        std::string_view config;
    };

#if defined(HASE_OPENPMD_PROBE_ADIOS)
    constexpr std::string_view ADIOS_CONFIG = R"({"backend":"adios2"})";
#endif
#if defined(HASE_OPENPMD_PROBE_HDF5)
    constexpr std::string_view HDF5_CONFIG = R"({"backend":"hdf5"})";
#endif
#if defined(HASE_OPENPMD_PROBE_ADIOS_SST)
    constexpr std::string_view SST_CONFIG
        = R"({"backend":"adios2","adios2":{"engine":{"type":"sst","parameters":{"DataTransport":"WAN","OpenTimeoutSecs":"1"}}}})";
#endif

    std::vector<Backend> candidates()
    {
        std::vector<Backend> result;
#if defined(HASE_OPENPMD_PROBE_ADIOS_SST)
        result.push_back({"adios-sst", ".sst", SST_CONFIG});
#endif
#if defined(HASE_OPENPMD_PROBE_ADIOS)
        result.push_back({"adios", ".bp", ADIOS_CONFIG});
#endif
#if defined(HASE_OPENPMD_PROBE_HDF5)
        result.push_back({"hdf5", ".h5", HDF5_CONFIG});
#endif
        return result;
    }

    void writeProbe(std::filesystem::path const& path, Backend const& backend)
    {
        io::Series series(path.string(), io::Access::CREATE_LINEAR, std::string{backend.config});
        auto iteration = series.writeIterations()[0];
        iteration.setAttribute("haseOpenPmdBackendProbe", std::string{backend.name});
        iteration.close();
        series.close();
    }

    bool readProbe(std::filesystem::path const& path, Backend const& backend)
    {
        io::Series series(path.string(), io::Access::READ_LINEAR, std::string{backend.config});
        for(auto iteration : series.readIterations())
        {
            auto value = iteration.getAttribute("haseOpenPmdBackendProbe").get<std::string>();
            series.close();
            return value == backend.name;
        }
        return false;
    }

    bool fileHandshake(Backend const& backend, std::filesystem::path const& path)
    {
        writeProbe(path, backend);
        return readProbe(path, backend);
    }

    bool streamingHandshake(Backend const& backend, std::filesystem::path const& path)
    {
        auto reader = std::async(std::launch::async, [&] { return readProbe(path, backend); });
        std::this_thread::sleep_for(std::chrono::milliseconds{100});
        auto writer = std::async(std::launch::async, [&] { writeProbe(path, backend); });
        writer.get();
        return reader.get();
    }

    bool handshake(Backend const& backend, std::filesystem::path const& dir)
    {
        auto path = dir / ("hase-openpmd-probe-" + std::string{backend.name} + std::string{backend.suffix});
        try
        {
            if(backend.name == "adios-sst")
            {
                return streamingHandshake(backend, path);
            }
            return fileHandshake(backend, path);
        }
        catch(std::exception const& exc)
        {
            std::cerr << "openPMD backend probe failed for " << backend.name << ": " << exc.what() << '\n';
            return false;
        }
    }
} // namespace

namespace hase::openpmd
{
    std::vector<std::string> availableOpenPmdBackends()
    {
        auto const workdir = std::filesystem::temp_directory_path()
                             / ("hase-openpmd-backend-probe-"
                                + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        std::filesystem::remove_all(workdir);
        std::filesystem::create_directories(workdir);

        std::vector<std::string> available;
        for(auto const& backend : candidates())
        {
            if(handshake(backend, workdir))
            {
                available.emplace_back(backend.name);
            }
        }
        std::filesystem::remove_all(workdir);
        return available;
    }

    void writeAvailableOpenPmdBackends(std::filesystem::path const& output)
    {
        std::filesystem::create_directories(output.parent_path());
        auto const available = availableOpenPmdBackends();
        if(available.empty())
        {
            throw std::runtime_error("No configured openPMD backend passed the create/read probe.");
        }

        std::ofstream stream(output);
        for(auto const& backend : available)
        {
            stream << backend << '\n';
        }
        if(!stream)
        {
            throw std::runtime_error("Failed to write openPMD backend probe output to " + output.string());
        }
    }
} // namespace hase::openpmd

#if defined(HASE_OPENPMD_BACKEND_PROBE_STANDALONE)
int main(int argc, char** argv)
{
    try
    {
        if(argc != 2)
        {
            std::cerr << "usage: " << argv[0] << " <output-file>\n";
            return 2;
        }
        hase::openpmd::writeAvailableOpenPmdBackends(std::filesystem::path{argv[1]});
        return 0;
    }
    catch(std::exception const& exc)
    {
        std::cerr << exc.what() << '\n';
        return 1;
    }
}
#endif
