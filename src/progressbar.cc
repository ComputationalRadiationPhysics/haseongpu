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

#include <cmath> /* log10 */
#include <iomanip> /* std::setfill, std::setw */
#include <fstream> /* std::ofstream */
#include <algorithm> /* std::max */
#include <chrono> /* steady_clock, time_point */
#include <sstream> /* stringstream */

#include <boost/filesystem/fstream.hpp> /* fs::fstream */
#include <boost/filesystem/path.hpp> /* fs::path */

#include <progressbar.hpp>
#include <logging.hpp>

namespace fs = boost::filesystem;
namespace chr = std::chrono;
void printWave(std::ostream& stream, unsigned tic, int progress, int length)
{
    for(int i = 0; i < progress; ++i)
    {
        switch((tic + i) % 12)
        {
            case 0:  stream << "ø"; break;
            case 1:  stream << "¤"; break;
            case 2:  stream << "º"; break;
            case 3:  stream << "°"; break;
            case 4:  stream << "`"; break;
            case 5:  stream << "°"; break;
            case 6:  stream << "º"; break;
            case 7:  stream << "¤"; break;
            case 8:  stream << "ø"; break;
            case 9:  stream << ","; break;
            case 10: stream << "¸"; break;
            case 11: stream << ","; break;
        }
    }

    for(int i = 0; i < length - progress; ++i)
    {
        stream << " ";
    }
}

void printBar(std::ostream& stream, int progress, int length)
{
    for(int i = 0; i < progress; ++i)
    {
        stream << "#";
    }
    for(int i = 0; i < length - progress; ++i)
    {
        stream << " ";
    }
}

std::string humanRepresentation(chr::duration<float> time)
{
    using seconds = chr::duration<int, std::ratio<1, 1>>;
    using minutes = chr::duration<int, std::ratio<60, 1>>;
    using hours   = chr::duration<int, std::ratio<3600, 1>>;

    auto const tSec  = chr::duration_cast<seconds>(time).count() % 60;
    auto const tMin  = chr::duration_cast<minutes>(time).count() % 60;
    auto const tHour = chr::duration_cast<hours>(time).count();

    std::stringstream ss;

    if(tHour) ss << tHour << "h ";
    if(tMin)  ss << tMin  << "m ";
    if(tSec)  ss << tSec  << "s";

    return ss.str();
}

ProgressBar::ProgressBar(){

}

void ProgressBar::reset()
{
    maxNTotal = 0;
    startTime = chr::steady_clock::now();
    part.store(0, std::memory_order_relaxed);   // atomic reset
    tic = 0;
    initialized = false;
}

void ProgressBar::printFancyProgressBar(unsigned nTotal)
{
#ifdef _WIN32
    const int length = 16;
#else
    const int length = 50;
#endif

    if(!initialized)
    {
        startTime = chr::steady_clock::now();
        initialized = true;
    }

    maxNTotal = std::max(maxNTotal, nTotal);

    unsigned const currentPart =
        part.fetch_add(1, std::memory_order_relaxed) + 1;

    if(maxNTotal == 0)
        return;

    auto const now = chr::steady_clock::now();
    chr::duration<float> const timeSpent = now - startTime;

    if(timeSpent.count() > 0.035f * tic || currentPart == maxNTotal)
    {
        ++tic;

        float const percentage =
            static_cast<float>(currentPart) / static_cast<float>(maxNTotal);

        auto const timeTotal =
            percentage > 0.0f
                ? timeSpent / percentage
                : chr::duration<float>::zero();

        auto const timeRemaining = timeTotal - timeSpent;

        unsigned const fillwidthPart =
            maxNTotal > 0
                ? static_cast<unsigned>(1 + std::log10(static_cast<double>(maxNTotal)))
                : 1;

        dout(V_PROGRESS | V_NOLABEL) << "\r";
        dout(V_PROGRESS) << "[";

#ifdef _WIN32
        printBar(
            dout(V_PROGRESS | V_NOLABEL),
            static_cast<int>(percentage * length),
            length);
#else
        printWave(
            dout(V_PROGRESS | V_NOLABEL),
            tic,
            static_cast<int>(percentage * length),
            length);
#endif

        dout(V_PROGRESS | V_NOLABEL) << "] ";
        dout(V_PROGRESS | V_NOLABEL)
            << std::setfill(' ') << std::setw(3)
            << static_cast<int>(percentage * 100) << "%";
        dout(V_PROGRESS | V_NOLABEL)
            << " (" << std::setfill(' ') << std::setw(fillwidthPart)
            << currentPart << "/" << maxNTotal << ")";
        dout(V_PROGRESS | V_NOLABEL)
            << " after " << humanRepresentation(timeSpent);
        dout(V_PROGRESS | V_NOLABEL)
            << " (" << humanRepresentation(timeTotal)
            << " total, " << humanRepresentation(timeRemaining)
            << " remaining)";
        dout(V_PROGRESS | V_NOLABEL) << std::flush;
    }
}
