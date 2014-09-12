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


#pragma once
#include <vector>

#define MAX_INTERPOLATION 1000
#define LAMBDA_START 905
#define LAMBDA_STOP 1095

std::vector<double> interpolateWavelength(const std::vector<double> sigma_y, const unsigned interpolation_range, const double lambda_start, const double lambda_stop);
std::vector<double> interpolateLinear(const std::vector<double> y, const std::vector<double> x, const unsigned nInterpolations);

