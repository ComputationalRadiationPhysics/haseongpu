#ifndef FILTERGRID_H
#define FILTERGRID_H 

#include "grid.h"

/**
 * @brief filter all geometric objects by whether they "might" be hit by the ray
 **/
__device__ void filter(const Grid *grid, const RayCu *ray, int* results);

#endif /* FILTERGRID_H */
