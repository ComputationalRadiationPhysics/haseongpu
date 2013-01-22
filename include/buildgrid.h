#ifndef BUILDGRID_H
#define BUILDGRID_H 

#include "grid.h"

/**
 * @brief calculate bounding boxes for numPrims prisms
 **/
void calcAabbs(const PrismCu *prisms, Aabb *aabbs, int numPrisms);

/**
 * @brief build grid based on the geometry-independent bounding boxes
 **/
void buildGrid(Grid *grid, const Aabb *aabbs, int numAabbs);

#endif /* BUILDGRID_H */
