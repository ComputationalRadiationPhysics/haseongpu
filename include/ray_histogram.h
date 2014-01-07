#ifndef RAY_HISTOGRAM_H
#define RAY_HISTOGRAM_H

#include <vector>

void ray_histogram(const std::vector<unsigned> totalRays, const unsigned max, const double mseThreshold, const std::vector<double> mseValues);

#endif
