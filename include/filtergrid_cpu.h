#ifndef FILTERGRID_CPU_H
#define FILTERGRID_CPU_H 

#include "grid.h"

/**
 * @brief filter all geometric objects by whether they "might" be hit by the ray
 **/
int* filterCpu(const Grid *grid, const RayCu *ray, int* result_size);

#endif /* FILTERGRID_CPU_H */
