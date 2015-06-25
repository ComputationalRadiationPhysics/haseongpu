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
#include <iostream>
#include <string>

// Idea from http://stackoverflow.com/questions/3342726/c-print-out-enum-value-as-text
// Answer by user SigTerm
// Switched char* to std::string
//
std::ostream& operator<<(std::ostream& out, const DeviceMode value){
    std::string s;
#define PROCESS_VAL(p) case(p): s = #p; break;
    switch(value){
        PROCESS_VAL(NO_DEVICE_MODE);
        PROCESS_VAL(GPU_DEVICE_MODE);
        PROCESS_VAL(CPU_DEVICE_MODE);
    }
#undef PROCESS_VAL
    return out << s;
}

std::ostream& operator<<(std::ostream& out, const ParallelMode value){
    std::string s;
#define PROCESS_VAL(p) case(p): s = #p; break;
    switch(value){
        PROCESS_VAL(NO_PARALLEL_MODE);
        PROCESS_VAL(THREADED_PARALLEL_MODE);
        PROCESS_VAL(MPI_PARALLEL_MODE);
        PROCESS_VAL(GRAYBAT_PARALLEL_MODE);
    }
#undef PROCESS_VAL
    return out << s;
}
