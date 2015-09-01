/**
 * Copyright 2015 Erik Zenker, Carlchristian Eckert
 *
 * This file is part of HASEonGPU
 *
 * HASEonGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HASEonGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HASEonGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */


#include "types.hpp"
#include <string>

const std::string DeviceMode::NONE  = "no_device_mode";
const std::string DeviceMode::GPU   = "gpu";
const std::string DeviceMode::CPU   = "cpu";

const std::string ParallelMode::NONE        = "no_parallel_mode";
const std::string ParallelMode::THREADED    = "threaded";

#if defined(MPI_FOUND)
const std::string ParallelMode::MPI         = "mpi";
#endif

#if defined(BOOST_MPI_FOUND) || defined(ZMQ_FOUND)
const std::string ParallelMode::GRAYBAT     = "graybat";
#endif

const std::string CompSwitch::parallel_mode     = "parallel-mode";
const std::string CompSwitch::device_mode       = "device-mode";
const std::string CompSwitch::ngpus             = "ngpus";
const std::string CompSwitch::repetitions       = "repetitions";
const std::string CompSwitch::adaptive_steps    = "adaptive-steps";
const std::string CompSwitch::min_sample_i      = "min-sample-i";
const std::string CompSwitch::max_sample_i      = "max-sample-i";

const std::string ExpSwitch::input_path     = "input-path";
const std::string ExpSwitch::output_path    = "output-path";
const std::string ExpSwitch::min_rays       = "min-rays";
const std::string ExpSwitch::max_rays       = "max-rays";
const std::string ExpSwitch::mse            = "mse-threshold";
const std::string ExpSwitch::reflection     = "reflection";
const std::string ExpSwitch::spectral       = "spectral-resolution";
