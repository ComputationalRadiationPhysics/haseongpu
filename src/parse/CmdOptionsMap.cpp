/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <core/logging.hpp>
#include <core/types.hpp>
#include <parse/CmdOptionsMap.hpp>

#include <iostream>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <utility>

namespace hase::parse
{

    std::string_view requireValue(int& i, int argc, char** argv, std::string_view opt)
    {
        if(i + 1 >= argc)
        {
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

    void CmdOptionsMap::add(std::string name, std::string description, std::string defaultValue)
    {
        Option opt;
        opt.name = std::move(name);
        opt.description = std::move(description);
        opt.defaultValue = std::move(defaultValue);
        opt.hasValue = false;

        options.emplace(opt.name, std::move(opt));
    }

    bool CmdOptionsMap::isExplicitlySet(std::string const& name) const
    {
        auto it = options.find(name);
        if(it == options.end())
        {
            return false;
        }
        return it->second.hasValue;
    }

    bool CmdOptionsMap::contains(std::string const& name) const
    {
        return options.find(name) != options.end();
    }

    void CmdOptionsMap::set(std::string const& name, std::string value)
    {
        auto it = options.find(name);
        if(it == options.end())
        {
            throw std::runtime_error("Unknown option: " + name);
        }
        it->second.set(std::move(value));
    }

    void CmdOptionsMap::printDescription() const
    {
        std::cout << "------------------------------------------" << std::endl;
        std::cout << "[HELP] HASEonGPU Argument Description:" << std::endl;

        for(auto const& opt : options | std::views::values)
        {
            std::cout << "--" << opt.name << '\n';
            std::cout << "  Description : " << opt.description << '\n';
            std::cout << "  Default     : " << (opt.defaultValue.empty() ? "<none>" : opt.defaultValue) << '\n';
        }

        std::cout << "------------------------------------------" << std::endl;
    }

    void CmdOptionsMap::checkRequired()
    {
        if(!isExplicitlySet(hase::core::ExpSwitch::input_path))
        {
            throw std::runtime_error("Missing required option: --input-path / -i");
        }

        if(!isExplicitlySet(hase::core::ExpSwitch::output_path))
        {
            options[hase::core::ExpSwitch::output_path].set(
                options[hase::core::ExpSwitch::input_path].as<std::string>());
        }
    }

    void CmdOptionsMap::construct()
    {
        add(hase::core::ExpSwitch::input_path, "The path to a folder that contains the input files", "");
        add(hase::core::ExpSwitch::output_path, "The path to a folder that contains the output files", "");

        add(hase::core::ExpSwitch::min_rays, "The minimal number of rays to use for each sample point", "100000");
        add(hase::core::ExpSwitch::max_rays, "The maximal number of rays to use for each sample point", "100000");
        add(hase::core::ExpSwitch::mse, "The MSE threshold used for adaptive/repetitive sampling", "0.1");
        add(hase::core::ExpSwitch::reflection, "Use reflections or not", "true");
        add(hase::core::ExpSwitch::spectral, "The number of samples used to interpolate spectral intensities", "");
        add(hase::core::ExpSwitch::monochromatic,
            "Use one constant sigmaA/sigmaE pair and skip wavelength interpolation",
            "false");

        add(hase::core::CompSwitch::parallel_mode, "Set the preferred way of parallelization (mpi, single)", "single");
        add(hase::core::CompSwitch::backend, "Set the device-backend to run the computation.", "gpu");
        add(hase::core::CompSwitch::numDevices, "The maximum number of GPUs to be used on a single node", "1");
        add(hase::core::CompSwitch::repetitions,
            "The number of repetitions to try before the number of rays is increased by adaptive sampling",
            "4");
        add(hase::core::CompSwitch::adaptive_steps,
            "The number of adaptive sampling steps used to split the range between min_rays and max_rays",
            "5");
        add(hase::core::CompSwitch::min_sample_i,
            "The minimum index of sample points to simulate",
            std::to_string(std::numeric_limits<unsigned>::max()));
        add(hase::core::CompSwitch::max_sample_i,
            "The maximal index of sample points to simulate",
            std::to_string(std::numeric_limits<unsigned>::max()));
        add(hase::core::CompSwitch::write_vtk, "Write VTK files of the computed ASE values", "false");
        add(hase::core::CompSwitch::backend, "The device-backend to use in order to perform the computation", "false");
        add("nPerNode",
            "Accepted for compatibility with wrapper scripts. Actual active ranks per node are detected at runtime.",
            "0");
        add("verbosity",
            "Set the verbosity levels",
            std::to_string(V_ERROR | V_INFO | V_WARNING | V_PROGRESS | V_STAT));
        add("config", "Location of an optional config file", "");
        add("help", "Print this help message and exit", "false");
    }

} // namespace hase::parse
