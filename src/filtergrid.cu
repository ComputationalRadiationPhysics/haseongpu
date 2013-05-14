#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "filtergrid.h"
#include "geometry_gpu.h"
#include "dynintset.h"

__device__ int cellIdx(const Grid *grid, int i, int j, int k){
	return k*grid->dim.x*grid->dim.y+j*grid->dim.x+i;
}

__device__ int iIdx(const Grid *grid, float x){
	//assert(x >= grid->aabb.min.x && x <= grid->aabb.max.x);
	
	return (x - grid->aabb.min.x)/(grid->aabb.max.x - grid->aabb.min.x) * grid->dim.x;
}

__device__ int jIdx(const Grid *grid, float y){
	//assert(y >= grid->aabb.min.y && y <= grid->aabb.max.y);
	
	return (y - grid->aabb.min.y)/(grid->aabb.max.y - grid->aabb.min.y) * grid->dim.y;
}

__device__ int kIdx(const Grid *grid, float z){
	//assert(z >= grid->aabb.min.z && z <= grid->aabb.max.z);
	
	return (z - grid->aabb.min.z)/(grid->aabb.max.z - grid->aabb.min.z) * grid->dim.z;
}

__device__ PointCu calcPointOnRay(const RayCu *ray, float t) {
  PointCu p;

  p.x = ray->P.x + t*ray->direction.x;
  p.y = ray->P.y + t*ray->direction.y;
  p.z = ray->P.z + t*ray->direction.z;
  
  return p;
}

__device__ int4 calcCellpos(const Grid *grid, PointCu p) {
  int4 result;

  result.x = ((int) p.x) % grid->dim.x;
  result.y = ((int) p.y) % grid->dim.y;
  result.z = ((int) p.z) % grid->dim.z;

  return result;
}

__device__ bool nextIntersection(const Grid *grid, const RayCu *ray, float *t, PointCu p) {
  float4 nextTs;
  nextTs.x = 0; nextTs.y = 0; nextTs.z = 0; // mark as invalid

  int step, next;
  float nextCoord, nextT;
  if(ray->direction.x != 0 && grid->dim.x != 0) {
    if(ray->direction.x > 0) step = 1;
    else step = -1;
    next = iIdx(grid, p.x) + step;
    nextCoord = ((float)next)/grid->dim.x * (grid->aabb.max.x - grid->aabb.min.x) + grid->aabb.min.x;
    nextT = (nextCoord - ray->P.x) / ray->direction.x;
    nextTs.x = nextT;
  }
  if(ray->direction.y != 0 && grid->dim.y != 0) {
    if(ray->direction.y > 0) step = 1;
    else step = -1;
    next = jIdx(grid, p.y) + step;
    nextCoord = ((float)next)/grid->dim.y * (grid->aabb.max.y - grid->aabb.min.y) + grid->aabb.min.y;
    nextT = (nextCoord - ray->P.y) / ray->direction.y;
    nextTs.y = nextT;
  }
  if(ray->direction.z != 0 && grid->dim.z != 0) {
    if(ray->direction.z > 0) step = 1;
    else step = -1;
    next = kIdx(grid, p.z) + step;
    nextCoord = ((float)next)/grid->dim.z * (grid->aabb.max.z - grid->aabb.min.z) + grid->aabb.min.z;
    nextT = (nextCoord - ray->P.z) / ray->direction.z;
    nextTs.z = nextT;
  }

  if(nextTs.x <= *t) nextTs.x = 1.1;
  if(nextTs.y <= *t) nextTs.y = 1.1;
  if(nextTs.z <= *t) nextTs.z = 1.1;
  *t = fmin(nextTs.x, fmin(nextTs.y, nextTs.z));
  if(*t > 1) return false;
  
  return true;
}

__device__ DynIntSet filter(const Grid *grid, const RayCu *ray) {
  DynIntSet results(256); // dynamically growing array with unique entries

  float t = 0;
  PointCu p;
  do {
    p = calcPointOnRay(ray, t);
    int4 cellpos = calcCellpos(grid, p);
    GridCell cell = grid->cellList[cellIdx(grid, cellpos.x, cellpos.y, cellpos.z)];
    results.push(cell.prismIdxList, cell.length);
  } while(nextIntersection(grid, ray, &t, p)); // changes t, but not p (is being done in next loop)

  return results;
}
