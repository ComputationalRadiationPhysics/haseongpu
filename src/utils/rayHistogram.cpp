/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 * Copyright 2026 Tim Hanel
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


#include <alpaka/alpaka.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

#ifdef _WIN32
#    include <Windows.h>
#endif

#include <core/logging.hpp>

namespace hase::utils
{

    void ray_histogram(
        std::vector<unsigned> const totalRays,
        unsigned const max,
        double const mseThreshold,
        std::vector<double> const mseValues,
        std::vector<unsigned> const droppedRays)
    {
        // length of the maximum number of samples (e.g. max==4210)
        int fillwidth = alpaka::math::log10(max) + 4;

        // maximum length of the filling bar
        unsigned maxLength = 50;

        // necessary size of the histogram
        std::map<unsigned, unsigned> histGreen;
        std::map<unsigned, unsigned> histRed;
        // if the entry doesn't exist, create it
        for(unsigned j = 0; j < totalRays.size(); ++j)
        {
            // std::map<unsigned,unsigned>::iterator it = hist.find(totalRays.at(j));
            if(histGreen.find(totalRays.at(j)) == histGreen.end())
            {
                histGreen.insert(std::pair<unsigned, unsigned>(totalRays.at(j), 0));
                histRed.insert(std::pair<unsigned, unsigned>(totalRays.at(j), 0));
            }
            if(mseValues.at(j) <= mseThreshold && droppedRays.at(j) == 0u)
            {
                histGreen.find(totalRays.at(j))->second++;
                // itG->second++;
            }
            else
            {
                histRed.find(totalRays.at(j))->second++;
                // itR->second++;
            }
        }


        std::map<unsigned, unsigned>::iterator itG;
        std::map<unsigned, unsigned>::iterator itR;
#ifndef _WIN32
        for(itG = histGreen.begin(), itR = histRed.begin(); itG != histGreen.end(); ++itG, ++itR)
        {
            hase::core::dout(V_STAT) << std::setw(fillwidth) << std::setfill(' ') << itG->first << " (";
            hase::core::dout(V_STAT | V_NOLABEL)
                << "\033[0;32m" << std::setw(alpaka::math::log10(totalRays.size()) + 3) << itG->second << "x";
            hase::core::dout(V_STAT | V_NOLABEL) << "\033[0m" << " / ";
            hase::core::dout(V_STAT | V_NOLABEL)
                << "\033[0;31m" << std::setw(alpaka::math::log10(totalRays.size()) + 3) << itR->second << "x";
            hase::core::dout(V_STAT | V_NOLABEL) << "\033[0m" << "):";

            // set color = green
            hase::core::dout(V_STAT | V_NOLABEL) << "\033[0;32m";
            for(unsigned j = 0; j < alpaka::math::ceil(maxLength * (float(itG->second) / totalRays.size())); ++j)
            {
                hase::core::dout(V_STAT | V_NOLABEL) << "#";
            }

            // set color = red
            hase::core::dout(V_STAT | V_NOLABEL) << "\033[0;31m";
            for(unsigned j = 0; j < alpaka::math::ceil(maxLength * (float(itR->second) / totalRays.size())); ++j)
            {
                hase::core::dout(V_STAT | V_NOLABEL) << "#";
            }
            hase::core::dout(V_STAT | V_NOLABEL) << std::endl;
        }

#else
        HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
        for(itG = histGreen.begin(), itR = histRed.begin(); itG != histGreen.end(); ++itG, ++itR)
        {
            hase::core::dout(V_STAT) << std::setw(fillwidth) << std::setfill(' ') << itG->first << " (";
            SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN | FOREGROUND_INTENSITY);
            hase::core::dout(V_STAT | V_NOLABEL)
                << std::setw(alpaka::math::log10(totalRays.size()) + 3) << itG->second << "x";
            SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN | FOREGROUND_RED | FOREGROUND_BLUE);
            hase::core::dout(V_STAT | V_NOLABEL) << " / ";
            SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_INTENSITY);
            hase::core::dout(V_STAT | V_NOLABEL)
                << std::setw(alpaka::math::log10(totalRays.size()) + 3) << itR->second << "x";
            SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN | FOREGROUND_RED | FOREGROUND_BLUE);
            hase::core::dout(V_STAT | V_NOLABEL) << "):";

            // set color = green
            SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN | FOREGROUND_INTENSITY);
            for(unsigned j = 0; j < alpaka::math::ceil(maxLength * (float(itG->second) / totalRays.size())); ++j)
            {
                hase::core::dout(V_STAT | V_NOLABEL) << "#";
            }

            // set color = red
            SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_INTENSITY);
            for(unsigned j = 0; j < alpaka::math::ceil(maxLength * (float(itR->second) / totalRays.size())); ++j)
            {
                hase::core::dout(V_STAT | V_NOLABEL) << "#";
            }
            hase::core::dout(V_STAT | V_NOLABEL) << std::endl;
        }
#endif
    }

} // namespace hase::utils
