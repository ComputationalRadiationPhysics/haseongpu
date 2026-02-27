/**
 * Copyright 2013-2026 Erik Zenker, Carlchristian Eckert, Marius Melzer, Tim Hanel
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


#include <vector>
#include <assert.h>

#include <interpolation.hpp>
#include <logging.hpp>

/**
 * @brief Returns the index of an value in vector v,
 *        that is smaller than t.
 *        
 * @param v vector
 * @param t bigger value
 * @return index of smaller value
 *         otherwise 0
 *
 */
unsigned getNextSmallerIndex(std::vector<double> v, double t){
  unsigned index = 0;
  for(unsigned i = 0; i < v.size(); ++i){
    if(v.at(i) < t) index = i;
    else break;
  }
  return index;
}

/**
 * @brief Returns the index of an value in vector v,
 *        that is bigger than t.
 *        
 * @param v vector
 * @param t smaller value
 * @return index of smaller value
 *         otherwise 0
 *
 */
unsigned getNextBiggerIndex(std::vector<double> v, double t){
  for(unsigned i = 0; i < v.size(); ++i){
    if(v.at(i) > t)
      return i;
  }
  return 0;
}

template <class T>
bool isStrictlySorted(const std::vector<T> &v){
    if(v.size() == 1) 
	return true;
    
    for(unsigned i = 0; i < v.size() - 1; ++i){
	if(v[i] >= v[i + 1]) {
	    return false;
	}

    }

    return true;
}


/**
 * @brief Interpolates the values of y to n values (nInterpolations) linear.
 *        The value of y[i] corresponds to the value of x[i].
 *        For Example could you interpolate 100 y values to 1000 y values,
 *        to reach a better resolution.
 *        The resulting x coordinates correspond to a equidistant spacing
 *
 * @param y                   y values
 * @param x                   x values of y (monoton increasing)
 * @param nInterpolations     number of interpolated values
 *
* @return A vector containing the linearly interpolated y values with size
 *         nInterpolations. The associated x-coordinates are distributed
 *         equidistantly over the interval (x[0],x[x.size()-1]).
 */
std::vector<double> interpolateLinear(const std::vector<double>& y, const std::vector<double>& x, const unsigned nInterpolations){
    assert(!x.empty());
    assert(x.size() == y.size());
    assert(nInterpolations >= y.size());
    assert(isStrictlySorted(x));
    if (nInterpolations == 0) return {};
    if (y.size() == 1) return std::vector<double>(nInterpolations, y.front());

    // Determine range of x values
    const double x_min = x.at(0);
    const double x_max = x.at(x.size() - 1);
    const double x_range = x_max - x_min;

    // Monochromatic case
    if(y.size() == 1){
	    return std::vector<double>(nInterpolations, y.front());
    }

    // Check x data if monoton increasing
    assert(isStrictlySorted(x));

    // Create equidistant interpolation on x axis
    std::vector<double> interpolated_x(nInterpolations, 0);
    for (unsigned i = 0; i < nInterpolations; ++i) {
      interpolated_x[i] = x_min + (static_cast<double>(i) * x_range) /
                                 static_cast<double>(nInterpolations - 1);
    }

    // Start to interpolate y values for every "new" x value
    std::vector<double> interpolated_y(nInterpolations, 0);
    unsigned j=0;
    for(unsigned i = 0; i < x.size()-1; ++i){
    	double x_orig=x[i];
    	double x_next=x[i+1];
    	double y_orig=y[i];
    	double y_next=y[i+1];
    	double m = (y_next - y_orig) / (x_next - x_orig);
    	double b = y_orig - (m * x_orig);
    	while (j<interpolated_y.size()&&interpolated_x[j]<x_next) {
    		interpolated_y[j]=b+m*interpolated_x[j];
    		++j;
    	}
    }
    interpolated_y.back() = y.back();
    return interpolated_y;
}
