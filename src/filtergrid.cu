#include "filtergrid.h"
#include "geometry_gpu.h"

__device__ int cellIdx(const Grid *grid, int i, int j, int k){
	return k*grid->dim.x*grid->dim.y+j*grid->dim.x+i;
}

//__global__ void filter(const Grid *grid, const RayCu *ray, char *results) {
//  int gx = threadIdx.x;
//  int gy = threadIdx.y;
//  int gz = threadIdx.z;
//
//  PointCu a, b, c, d;
//  PrismCu p1, p2;
//
//  VectorCu gridAabbDim = createVector(grid->aabb.min, grid->aabb.max);
//  VectorCu cellDim;
//  cellDim.x = gridAabbDim.x / grid->dim.x;
//  cellDim.y = gridAabbDim.y / grid->dim.y;
//  cellDim.z = gridAabbDim.z / grid->dim.z;
//
//  a = createPoint(gx*cellDim.x, gy*cellDim.y, gz*cellDim.z, cellDim.z);
//  b = createPoint(gx*cellDim.x, gy*cellDim.y, gz*cellDim.z, cellDim.z);
//  c = createPoint(gx*cellDim.x, gy*cellDim.y, gz*cellDim.z, cellDim.z);
//  d = createPoint(gx*cellDim.x, gy*cellDim.y, gz*cellDim.z, cellDim.z);
//
//  p1.t1 = createTriangle(a, b, c);
//  p2.t1 = createTriangle(b, c, d);
//
//  results[cellIdx(grid, gx, gy, gz)] = (collide_prism_gpu(p1, *ray) != 0) || (collide_prism_gpu(p2, *ray) != 0);
//}

//int iIdx(const Grid *grid, float x){
//	assert(x >= grid->aabb.min.x && x <= grid->aabb.max.x);
//	
//	return (x - grid->aabb.min.x)/(grid->aabb.max.x - grid->aabb.min.x) * grid->dim.x;
//}
//
//int jIdx(const Grid *grid, float y){
//	assert(y >= grid->aabb.min.y && y <= grid->aabb.max.y);
//	
//	return (y - grid->aabb.min.y)/(grid->aabb.max.y - grid->aabb.min.y) * grid->dim.y;
//}
//
//int kIdx(const Grid *grid, float z){
//	assert(z >= grid->aabb.min.z && z <= grid->aabb.max.z);
//	
//	return (z - grid->aabb.min.z)/(grid->aabb.max.z - grid->aabb.min.z) * grid->dim.z;
//}
//
//void cellIdx2ijk(const Grid *grid, int cellIdx, int *i, int *j, int *k){
//	
//	int dimXdimY=grid->dim.x*grid->dim.y;
//	int tmp;
//	
//	*k = cellIdx / dimXdimY;
//	
//	tmp=cellIdx - *k * dimXdimY;
//	//*j = (cellIdx - *k * dimXdimY) / grid->dim.x;
//	*j = tmp / grid->dim.x;
//	
//	//*i = cellIdx - *k * dimXdimY - *j* grid->dim.x;
//	*i = tmp - *j * grid->dim.x;
//}

__device__ PointCu* nearest(PointCu p, PointCu *a, PointCu *b) {
  if(b->w == -1 || distance_gpu(p, *a) < distance_gpu(p, *b)) return a;
  return b;
}

__device__ int round(int number, int cellsize, float dir) {
  if(number % cellsize == 0) {
    if(dir > 0) return (number+cellsize);
    else if(dir < 0) return (number-cellsize);
    else return number;
  } else {
    int mod = number % cellsize;
    if(dir > 0) return number + (cellsize - mod);
    else return number - (cellsize + mod);
  }
}

__device__ bool nextIntersection(const RayCu *ray, PointCu p, int4 dim, PointCu *result) {
  PointCu end = addVectorToPoint(ray->P, ray->direction);
  if(end.x == p.x && end.y == p.y && end.z == p.z) return false;
  PointCu x, y, z;

  x.w = 0; y.w = 0; z.w = 0; // mark as invalid

  x.x = round(p.x);
  if(x.x != p.x) x.w = 1;
  y.y = round(p.y);
  if(x.x != p.x) y.w = 1;
  x.x = round(p.z);
  if(z.z != p.z) z.w = 1;

  result = nearest(p, nearest(p, &x, &y), &z);
  if(nearest(p, result, &end) == &end) return false; // check if result is more far than end
  return true;
}

__device__ int4 calcCellpos(const Grid *grid, PointCu p) {
  int4 result;

  result.x = ((int) p.x) % grid->dim.x;
  result.y = ((int) p.y) % grid->dim.y;
  result.z = ((int) p.z) % grid->dim.z;

  return result;
}

__device__ int* cuRealloc(int *old, int *size) {
  int oldSize = *size;
  int* res = (int*) malloc(oldSize * 2);
  *size *= 2;

  for(int i=0; i<oldSize; ++i)
    res[i] = old[i];

  free(old);
  return res;
}

__device__ void addPrismsToResult(const GridCell *cell, int* results, int *results_count, int *results_size) {
  int limit = *results_count + cell->length;
  if(limit >= *results_size) {
    while(limit >= *results_size) *results_size += 256;
    results = cuRealloc(results, results_size);
  }

  for(int i=0; i<cell->length; ++i)
    results[*results_count + i] = cell->prismIdxList[i];

  *results_count += cell->length;
}

__device__ int filter(const Grid *grid, const RayCu *ray, int* results) {
  results = (int*) malloc(256 * sizeof(int));
  int results_count = 0;
  int results_size = 256;

  PointCu p = ray->P;
  PointCu res;
  while(nextIntersection(ray, p, grid->dim, &res)) {
    int4 cellpos = calcCellpos(grid, p);
    GridCell cell = grid->cellList[cellIdx(grid, cellpos.x, cellpos.y, cellpos.z)];
    addPrismsToResult(&cell, results, &results_count, &results_size); // uses sometimes realloc and changes results, results_count and results_size
    p = res;
  }

  return results_count;
}
