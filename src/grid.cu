#include <cstdio>
#include <cassert>

#include "grid.h"

//???????????
//dynamic list???
//save??
int filter(const Grid *grid, const RayCu *ray, int* results);

//BoundingBox for a Prism
void calcAabbs(const PrismCu *prisms, Aabb *aabbs, int numPrisms);

void calcGridAabb(const Grid *grid, const Aabb *objectAabbs, int numAabbs);

//define Dimensions
void setDimGrid(Grid *grid, int4 dimension);
//or grid.dimX=...

//build Grid
//numCells
void buildGrid(Grid *grid, const Aabb *aabbs, int numAabbs);

//
void fillGrid(const Grid *grid, const Aabb *aabbs, int numAabbs);

//calcAabbs

//buildGrid
	//calcGridAabb
	//define Dimensions
  //fillGrid

int cellIdx(const Grid *grid, int i, int j, int k){
	return k*grid->dimX*grid->dimY+j*grid->dimX+i;
}

int iIdx(const Grid *grid, float x){
	assert(x >= grid->aabb.min.x && x <= grid->aabb.max.x);
	
	return (x - grid->aabb.min.x)/(grid->aabb.max.x - grid->aabb.min.x) * grid->dimX;
}

int jIdx(const Grid *grid, float y){
	assert(y >= grid->aabb.min.y && y <= grid->aabb.max.y);
	
	return (y - grid->aabb.min.y)/(grid->aabb.max.y - grid->aabb.min.y) * grid->dimY;
}

int kIdx(const Grid *grid, float z){
	assert(z >= grid->aabb.min.z && z <= grid->aabb.max.z);
	
	return (z - grid->aabb.min.z)/(grid->aabb.max.z - grid->aabb.min.z) * grid->dimZ;
}

void cellIdx2ijk(const Grid *grid, int cellIdx, int *i, int *j, int *k){
	
	int dimXdimY=grid->dimX*grid->dimY;
	int tmp;
	
	*k = cellIdx / dimXdimY;
	
	tmp=cellIdx - *k * dimXdimY;
	//*j = (cellIdx - *k * dimXdimY) / grid->dimX;
	*j = tmp / grid->dimX;
	
	//*i = cellIdx - *k * dimXdimY - *j* grid->dimX;
	*i = tmp - *j * grid->dimX;
}
