/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 *
 * This file is part of HASENonGPU
 *
 * HASENonGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HASENonGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HASENonGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */


/**
 * @author Carlchristian Eckert
 * @license GPLv3
 *
 * Library to define flexible verbosity based on a global variable
 */
#ifndef OCTRACE_LOGGING_H
#define OCTRACE_LOGGING_H

#include <iostream>

#define V_QUIET 0
#define V_ERROR 1
#define V_WARNING 2
#define V_INFO 4
#define V_STAT 8
#define V_DEBUG 16
#define V_PROGRESS 32
#define V_NOLABEL 64

#define COLOR_ERROR 31 //red
#define COLOR_WARN 31  //red
#define COLOR_INFO 32  //green
#define COLOR_STATISTIC 32 //green
#define COLOR_DEBUG 36  //cyan


extern unsigned verbosity;

struct nullstream : std::ostream {
  nullstream() : std::ostream(0) { }
};

/**
 * @brief An output stream that writes to std::cout or std::cerr (depending
 *        on the verbosity that is set.
 *
 *        usage: 
 *        dout(V_ERROR) << ... 
 *        writes only, if the verbosity-flag V_ERROR is activated
 *
 * @param activation_level a bitmask containing all the bits, on which 
 *        verbosity levels the output should appear
 *
 */
std::ostream& dout(unsigned activation_level);


#endif
