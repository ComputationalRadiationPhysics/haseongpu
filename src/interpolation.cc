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
 *
 * @param y                   y values
 * @param x                   x values of y (monoton increasing)
 * @param nInterpolations     number of interpolated values
 *
 * @return vector of linear interpolated y values with size of nInterpolations
 */
std::vector<double> interpolateLinear(const std::vector<double> y, const std::vector<double> x, const unsigned nInterpolations){
    assert(nInterpolations >= y.size());
    assert(y.size() == x.size());

    // Determine range of x values
    const double x_min = x.at(0);
    const double x_max = x.at(x.size() - 1);
    const double x_range = x_max - x_min;
    assert(y.size() >= x_range);

    // Monochromatic case
    if(y.size() == 1){
	return std::vector<double>(nInterpolations, y.at(0));
    }

    // Check x data if monoton increasing
    assert(isStrictlySorted(x));

    // Create equidistant interpolation on x axis
    std::vector<double> interpolated_x(nInterpolations, 0);
    for(unsigned i = 0; i < nInterpolations; ++i){
	interpolated_x.at(i) = x_min + (i * (x_range / nInterpolations));
    }

    // Start to interpolate y values for every x value
    std::vector<double> interpolated_y(nInterpolations, 0);
    for(unsigned i = 0; i < interpolated_x.size(); ++i){
	// Get index of points before and after x
	double y1_i = getNextSmallerIndex(x, interpolated_x.at(i));
	double y2_i = getNextBiggerIndex(x, interpolated_x.at(i));
	int y_diff = y2_i - y1_i;

	if(y_diff == 1){
	    // First point p1=(x1/y1) before x
	    double x1 = x_min + y1_i;
	    double y1 = y.at(y1_i);

	    // Second point p2=(x2/y2) after x
	    double x2 = x_min + y2_i;
	    double y2 = y.at(y2_i);
	    assert(y.size() >= y1_i);

	    // linear function between p1 and p2 (y=mx+b)
	    double m = (y2 - y1) / (x2 / x1);
	    double b = y1 - (m * x1);

	    // Interpolate y from linear function
	    interpolated_y.at(i) = m * interpolated_x.at(i) + b;

	}
	else if(y_diff == 2){
	    // No interpolation needed
	    interpolated_y.at(i) = y.at(y1_i + 1);
	}
	else {
	    dout(V_ERROR) << "Index of smaller and bigger sigma too seperated" << std::endl;
	    exit(0);
	}
    
    }
  
    return interpolated_y;
}



/************************************************************************************/
/************************************************************************************/

/**
 * @brief Interpolates the values of sigma_y to n values(interpolation range) linear. 
 *        With the assumption, they are distributed between lambda_start
 *        and lambda_stop equidistant. For Example could you interpolate
 *        100 sigma values to 1000 sigma values, to reach a better resolution.
 *
 * @deprecated use interpolateLinear instead !
 *
 * @param sigma_y             y values
 * @param interpolation_range number of interpolated values
 * @param lambda_start        start of x range
 * @param lambda_stop         stop of x range
 *
 * @return vector of linear interpolated values
 */
std::vector<double> interpolateWavelength(const std::vector<double> sigma_y, const unsigned interpolation_range, const double lambda_start, const double lambda_stop){
  assert(interpolation_range >= sigma_y.size());
  assert(lambda_stop >= lambda_start);

  // Monochromatic case
  if(sigma_y.size() == 1){
    return std::vector<double>(1, sigma_y.at(0));
  }

  std::vector<double> y(interpolation_range, 0);
  const double lambda_range = lambda_stop - lambda_start;
  assert(sigma_y.size() >= lambda_range);

  // Generate sigma_x
  std::vector<double> sigma_x;
  for(unsigned i = lambda_start; i <= lambda_stop; ++i){
    sigma_x.push_back(i);
  }
  
  for(unsigned i = 0; i < interpolation_range; ++i){
    double x = lambda_start + (i * (lambda_range / interpolation_range));

    // Get index of points before and after x
    double y1_i = getNextSmallerIndex(sigma_x, x);
    double y2_i = getNextBiggerIndex(sigma_x, x);
    int sigma_diff = y2_i - y1_i;

    if(sigma_diff == 1){
      // First point p1=(x1/y1) before x
      double x1 = lambda_start + y1_i;
      double y1 = sigma_y.at(y1_i);

      // Second point p2=(x2/y2) after x
      double x2 = lambda_start + y2_i;
      double y2 = sigma_y.at(y2_i);
      assert(sigma_y.size() >= y1_i);

      // linear function between p1 and p2 (y=mx+b)
      double m = (y2 - y1) / (x2 / x1);
      double b = y1 - (m * x1);

      // Interpolate y from linear function
      y.at(i) = m * x + b;

    }
    else if(sigma_diff == 2){
      // No interpolation needed
      y.at(i) = sigma_y.at(y1_i + 1);
    }
    else {
      dout(V_ERROR) << "Index of smaller and bigger sigma too seperated" << std::endl;
      exit(0);
    }
    
  }
  
  return y;
}
