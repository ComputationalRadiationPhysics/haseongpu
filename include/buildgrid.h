#ifndef BUILDGRID_H
#define BUILDGRID_H 

#include "grid.h"

/**
 * @brief calculate bounding boxes for numPrims prisms
 **/
void calcAabbs(const PrismCu *prisms, Aabb *aabbs, int numPrisms);

/**
 * @brief set dimension of the grid
 **/
void setDimGrid(Grid *grid, int numTrianglesPerLayer, int numLayers, int numPrismPerCell);

/**
 * @brief build grid based on the geometry-independent bounding boxes
 **/
void buildGrid(Grid *grid, const Aabb *aabbs, int numAabbs);

/**
 * @brief copy Grid to GPU
 **/
void copyGrid2GPU(const Grid *grid, Grid *d_grid);


/**
 * @brief copy Grid to GPU
 * d_grid is a host variable -> call by reference 
 * copy it to a kernel call by value after using this function
 **/
void prepareGrid(Grid *d_grid, int4 dimGrid, const PrismCu *prisms, int numPrisms);

#endif /* BUILDGRID_H */
