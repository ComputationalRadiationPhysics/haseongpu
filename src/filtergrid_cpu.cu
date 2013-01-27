#include "filtergrid.h"
#include "geometry.h"
#include "print.h"
#include <stdio.h>
#include "assert.h"
#include "unistd.h"
#include "math.h"

int cellIdxCpu(const Grid *grid, int i, int j, int k){
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

int iIdxCpu(const Grid *grid, float x){
	//assert(x >= grid->aabb.min.x && x <= grid->aabb.max.x);
	
	return (x - grid->aabb.min.x)/(grid->aabb.max.x - grid->aabb.min.x) * grid->dim.x;
}

int jIdxCpu(const Grid *grid, float y){
	//assert(y >= grid->aabb.min.y && y <= grid->aabb.max.y);
	
	return (y - grid->aabb.min.y)/(grid->aabb.max.y - grid->aabb.min.y) * grid->dim.y;
}

int kIdxCpu(const Grid *grid, float z){
	//assert(z >= grid->aabb.min.z && z <= grid->aabb.max.z);
	
	return (z - grid->aabb.min.z)/(grid->aabb.max.z - grid->aabb.min.z) * grid->dim.z;
}

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

PointCu calcPointOnRayCpu(const RayCu *ray, float t) {
  printf("calcPointOnRay: %f: \n", t);
  PointCu p;

  p.x = ray->P.x + t*ray->direction.x;
  p.y = ray->P.y + t*ray->direction.y;
  p.z = ray->P.z + t*ray->direction.z;
  
  return p;
}

bool nextIntersectionCpu(const Grid *grid, const RayCu *ray, float *t, PointCu p) {
  puts("In nextIntersectionCpu");

  float4 nextTs;
  nextTs.x = 0; nextTs.y = 0; nextTs.z = 0; // mark as invalid

  int step, next;
  float nextCoord, nextT;
  if(ray->direction.x != 0 && grid->dim.x != 0) {
    if(ray->direction.x > 0) step = 1;
    else step = -1;
    next = iIdxCpu(grid, p.x) + step;
    nextCoord = ((float)next)/grid->dim.x * (grid->aabb.max.x - grid->aabb.min.x) + grid->aabb.min.x;
    nextT = (nextCoord - ray->P.x) / ray->direction.x;
    printf("in x: step: %i, next: %i, nextCoord: %f, nextT: %f\n", step, next, nextCoord, nextT);
    nextTs.x = nextT;
  }
  if(ray->direction.y != 0 && grid->dim.y != 0) {
    if(ray->direction.y > 0) step = 1;
    else step = -1;
    next = jIdxCpu(grid, p.y) + step;
    nextCoord = ((float)next)/grid->dim.y * (grid->aabb.max.y - grid->aabb.min.y) + grid->aabb.min.y;
    nextT = (nextCoord - ray->P.y) / ray->direction.y;
    printf("in y: step: %i, next: %i, nextCoord: %f, nextT: %f\n", step, next, nextCoord, nextT);
    nextTs.y = nextT;
  }
  if(ray->direction.z != 0 && grid->dim.z != 0) {
    if(ray->direction.z > 0) step = 1;
    else step = -1;
    next = kIdxCpu(grid, p.z) + step;
    nextCoord = ((float)next)/grid->dim.z * (grid->aabb.max.z - grid->aabb.min.z) + grid->aabb.min.z;
    nextT = (nextCoord - ray->P.z) / ray->direction.z;
    printf("in z: step: %i, next: %i, nextCoord: %f, nextT: %f\n", step, next, nextCoord, nextT);
    nextTs.z = nextT;
  }

  if(nextTs.x <= *t) nextTs.x = 1.1;
  if(nextTs.y <= *t) nextTs.y = 1.1;
  if(nextTs.z <= *t) nextTs.z = 1.1;
  *t = fmin(nextTs.x, fmin(nextTs.y, nextTs.z));
  if(*t > 1) return false;
  
  return true;
}

int4 calcCellposCpu(const Grid *grid, PointCu p) {
  puts("In calcCellposCpu");
  int4 result;

  result.x = ((int) p.x) % grid->dim.x;
  result.y = ((int) p.y) % grid->dim.y;
  result.z = ((int) p.z) % grid->dim.z;

  return result;
}

int* cuReallocCpu(int *old, int *size) {
  puts("In cuReallocCpu");
  int oldSize = *size;
  int* res = (int*) malloc(oldSize * 2);
  *size *= 2;

  for(int i=0; i<oldSize; ++i)
    res[i] = old[i];

  free(old);
  return res;
}

bool isDuplicateCpu(int *a, int length, int o) {
  for(int i=0; i<length; ++i) {
    if(a[i] == o) return true;
  }
  return false;
}

void addPrismsToResultCpu(const GridCell *cell, int* results, int *results_count, int *results_size) {
  puts("In addPrismsToResultCpu");
  int limit = *results_count + cell->length;
  if(limit >= *results_size) {
    while(limit >= *results_size) *results_size += 256;
    results = cuReallocCpu(results, results_size);
  }

  int added = 0;
  for(int i=0; i<cell->length; ++i) {
    int current = cell->prismIdxList[i];
    if(!isDuplicateCpu(results, *results_count, current)) {
      printf("save %i to %i\n", current, *results_count + i);
      results[*results_count + added] = current;
      added++;
    }
  }
  *results_count += added;
}

int* filterCpu(const Grid *grid, const RayCu *ray, int *result_size) {
  puts("In filterCpu");
  int *results = (int*) malloc(256 * sizeof(int));
  int results_count = 0;
  int results_size = 256;

  float t = 0;
  PointCu p;
  do {
    puts("");
    sleep(1);
    p = calcPointOnRayCpu(ray, t);
    int4 cellpos = calcCellposCpu(grid, p);
    GridCell cell = grid->cellList[cellIdxCpu(grid, cellpos.x, cellpos.y, cellpos.z)];
    addPrismsToResultCpu(&cell, results, &results_count, &results_size); // uses sometimes realloc and changes results, results_count and results_size
  } while(nextIntersectionCpu(grid, ray, &t, p)); // changes t, but not p (is being done in next loop)

  *result_size = results_count;
  return results;
}
