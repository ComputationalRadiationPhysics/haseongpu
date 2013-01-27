#ifndef FILTERGRID_H
#define FILTERGRID_H 

#include "grid.h"

/**
 * @brief filter all geometric objects by whether they "might" be hit by the ray
 **/
__device__ int* filter(const Grid *grid, const RayCu *ray, int* result_size);

#endif /* FILTERGRID_H */
