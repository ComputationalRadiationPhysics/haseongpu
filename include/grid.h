#ifndef GRID_H
#define GRID_H 

#include "datatypes.h"

//***************
// DATASTRUCTURES 
//***************

//Axis Aligned Bounding Box
struct Aabb {
	PointCu min;
	PointCu max;
};

struct GridCell {
	int* prismIdxList;
	int length;
};

struct Grid {
	Aabb aabb;
	GridCell* cellList;
  int4 dimensions;
};

//***********
// BUILD GRID
//***********

/**
 * @brief calculate bounding boxes for numPrims prisms
 **/
void calcAabbs(const PrismCu *prisms, Aabb *aabbs, int numPrisms);

/**
 * @brief build grid based on the geometry-independent bounding boxes
 **/
void buildGrid(Grid *grid, const Aabb *aabbs, int numAabbs);

//************
// SEARCH GRID
//************

/**
 * @brief filter all geometric objects by whether they "might" be hit by the ray
 **/
int filter(const Grid *grid, const RayCu *ray, int* results);

#endif /* GRID_H */
