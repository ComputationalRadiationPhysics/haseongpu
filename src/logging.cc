/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
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

#include <iostream>

#ifdef _WIN32
#include <Windows.h>
#endif

#include <logging.hpp>


std::ostream& dout(unsigned activation_level) {
  static nullstream dummy;

  if(!(activation_level & verbosity)){
    return dummy;
  }

  if(activation_level & V_NOLABEL){
    return std::cout;
  }

  if(activation_level & V_ERROR){
#ifdef _WIN32
    HANDLE hConsole = GetStdHandle (STD_OUTPUT_HANDLE);
    SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_INTENSITY);
    return std::cerr << "[ERROR] ";
#else
    return std::cerr << "\033[0;" << COLOR_ERROR << "m[ERROR] ";
#endif
  }

  if(activation_level & V_WARNING){
#ifdef _WIN32
    HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_INTENSITY);
    return std::cerr << "[WARNING] ";
#else
    return std::cerr << "\033[0;" << COLOR_WARN << "m[WARNING] ";
#endif
  }

  if(activation_level & V_INFO){
#ifdef _WIN32
    HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleTextAttribute(hConsole, FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED);
    return std::cerr << "[INFO] ";
#else
    //return std::cout << "\033[0;" << COLOR_INFO << "m[INFO] ";
    return std::cout << "\033[0" << "m[INFO] ";
#endif
  }

  if(activation_level & V_STAT){
#ifdef _WIN32
    HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleTextAttribute(hConsole, FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED);
    return std::cerr << "[STATISTIC] ";
#else
    //return std::cout << "\033[0;" << COLOR_STATISTIC << "m[STATISTIC] ";
    return std::cout << "\033[0" << "m[STATISTIC] ";
#endif
  }

  if(activation_level & V_PROGRESS){
#ifdef _WIN32
    HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleTextAttribute(hConsole, FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED);
    return std::cerr << "[PROGRESS] ";
#else
    return std::cout << "\033[0" << "m[PROGRESS] ";
#endif
  }

  if(activation_level & V_DEBUG){
#ifdef _WIN32
    HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleTextAttribute(hConsole, FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_INTENSITY);
    return std::cerr << "[DEBUG] ";
#else
    return std::cerr << "\033[0;" << COLOR_DEBUG << "m[DEBUG] ";
#endif
  }


  return std::cout;
}

