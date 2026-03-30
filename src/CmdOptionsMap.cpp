/**
* Copyright 2026 Tim Hanel
*/
#include "CmdOptionsMap.hpp"

#include <iostream>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <utility>

std::string_view requireValue(int& i, int argc, char** argv, std::string_view opt)
{
    if (i + 1 >= argc) {
        throw std::runtime_error("Missing value for option " + std::string(opt));
    }
    ++i;
    return argv[i];
}

void CmdOptionsMap::Option::set(std::string newValue)
{
    value = std::move(newValue);
    hasValue = true;
}

CmdOptionsMap::CmdOptionsMap()
{
    construct();
}

void CmdOptionsMap::add(
    std::string name,
    std::string description,
    std::string defaultValue)
{
    Option opt;
    opt.name = std::move(name);
    opt.description = std::move(description);
    opt.defaultValue = std::move(defaultValue);
    opt.hasValue = false;

    options.emplace(opt.name, std::move(opt));
}

bool CmdOptionsMap::isExplicitlySet(const std::string& name) const
{
    auto it = options.find(name);
    if (it == options.end()) {
        return false;
    }
    return it->second.hasValue;
}

bool CmdOptionsMap::contains(const std::string& name) const
{
    return options.find(name) != options.end();
}

void CmdOptionsMap::set(const std::string& name, std::string value)
{
    auto it = options.find(name);
    if (it == options.end()) {
        throw std::runtime_error("Unknown option: " + name);
    }
    it->second.set(std::move(value));
}

void CmdOptionsMap::printDescription() const
{
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "[HELP] HASEonGPU Argument Description:" << std::endl;

    for (const auto& opt : options | std::views::values) {
        std::cout << "--" << opt.name << '\n';
        std::cout << "  Description : " << opt.description << '\n';
        std::cout << "  Default     : "
                  << (opt.defaultValue.empty() ? "<none>" : opt.defaultValue)
                  << '\n';
    }

    std::cout << "------------------------------------------" << std::endl;
}

void CmdOptionsMap::checkRequired()
{
    if (!isExplicitlySet(ExpSwitch::input_path)) {
        throw std::runtime_error("Missing required option: --input-path / -i");
    }

    if (!isExplicitlySet(ExpSwitch::output_path)) {
        options[ExpSwitch::output_path].set(options[ExpSwitch::input_path].as<std::string>());
    }
}

void CmdOptionsMap::construct()
{
    add(ExpSwitch::input_path,    "The path to a folder that contains the input files", "");
    add(ExpSwitch::output_path,   "The path to a folder that contains the output files", "");

    add(ExpSwitch::min_rays,      "The minimal number of rays to use for each sample point", "100000");
    add(ExpSwitch::max_rays,      "The maximal number of rays to use for each sample point", "100000");
    add(ExpSwitch::mse,           "The MSE threshold used for adaptive/repetitive sampling", "0.1");
    add(ExpSwitch::reflection,    "Use reflections or not", "true");
    add(ExpSwitch::spectral,      "The number of samples used to interpolate spectral intensities", "");

    add(CompSwitch::parallel_mode, "Set the preferred way of parallelization (mpi, graybat, threaded), only valid with --device_mode=gpu", "threaded");
    add(CompSwitch::device_mode,   "Set the device to run the calculation (cpu, gpu)", "gpu");
    add(CompSwitch::ngpus,         "The maximum number of GPUs to be used on a single node", "1");
    add(CompSwitch::repetitions,   "The number of repetitions to try before the number of rays is increased by adaptive sampling", "4");
    add(CompSwitch::adaptive_steps,"The number of adaptive sampling steps used to split the range between min_rays and max_rays", "5");
    add(CompSwitch::min_sample_i,  "The minimum index of sample points to simulate", std::to_string(std::numeric_limits<unsigned>::max()));
    add(CompSwitch::max_sample_i,  "The maximal index of sample points to simulate", std::to_string(std::numeric_limits<unsigned>::max()));
    add(CompSwitch::write_vtk,     "Write VTK files of the computed ASE values", "false");

    add("verbosity",               "Set the verbosity levels", std::to_string(V_ERROR | V_INFO | V_WARNING | V_PROGRESS | V_STAT));
    add("config",                  "Location of an optional config file", "");
    add("help",                    "Print this help message and exit", "false");
}