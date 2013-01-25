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
  int4 dim;
};

#endif /* GRID_H */
