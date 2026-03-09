/**
 * Copyright 2015 Erik Zenker, Carlchristian Eckert, Marius Melzer
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


/**
 * @brief A collection of functions to print the progress
 *        of a process as a nice progressbar
 *
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @licence GPLv3
 *
 **/

#pragma once

#pragma once

#include <atomic>
#include <chrono>
#include <iosfwd>
#include <string>
namespace chr = std::chrono;

struct ProgressBar
{
    unsigned maxNTotal = 0;

    chr::time_point<chr::steady_clock> startTime{};
    std::atomic<unsigned> part{0};
    unsigned tic = 0;
    bool initialized = false;

    ProgressBar();

    void reset();
    void printFancyProgressBar(unsigned nTotal);
};
