/**
 * Copyright 2014 Erik Zenker, Carlchristian Eckert, Marius Melzer
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
 * @author Carlchristian Eckert
 * @licence GPLv3
 *
 */

#pragma once

/**
 * @brief this allows the use of isnan() for int and unsigned
 * in the template function fileToVector()
 */
inline bool isNaN(const int){
  return false;
}

inline bool isNaN(const unsigned int){
  return false;
}

inline bool isNaN(const float i){
  return isnan(i);
}

inline bool isNaN(const double i){
  return isnan(i);
}

inline bool isNaN(const long double i){
  return isnan(i);
}

