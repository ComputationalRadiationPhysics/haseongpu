#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <vector>

#define MAX_INTERPOLATION 1000
#define LAMBDA_START 905
#define LAMBDA_STOP 1095

std::vector<double> interpolateWavelength(const std::vector<double> sigma_y, const unsigned interpolation_range, const double lambda_start, const double lambda_stop);

#endif /* INTERPOLATION_H */
