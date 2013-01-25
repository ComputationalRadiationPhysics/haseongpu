#include "grid.h"

//allocate memory for aabbs before start the function
//assumption (Annahme): triangles in a plane with z=const
void calcAabbs(const PrismCu *prisms, Aabb *aabbs, int numPrisms)
{
  int i;
  float xmin,xmax,ymin,ymax;

  Aabb *pAabb;
  const TriangleCu *ptriangle;

  for(i=0;i<numPrisms;++i){
    //set trianglePointer
    ptriangle=&prisms[i].t1;

    //start values for xmin,xmax,ymin,ymax from first Point of the triangle
    xmax=xmin=ptriangle->A.x;
    ymax=ymin=ptriangle->A.y;

    if( ptriangle->B.x > xmax) xmax = ptriangle->B.x;
    else if( ptriangle->B.x < xmin) xmin = ptriangle->B.x;

    if( ptriangle->C.x > xmax) xmax = ptriangle->C.x;
    else if( ptriangle->C.x < xmin) xmin = ptriangle->C.x;

    if( ptriangle->B.y > ymax) ymax = ptriangle->B.y;
    else if( ptriangle->B.y < ymin) ymin = prisms[i].t1.B.y;

    if( ptriangle->C.y > ymax) ymax = ptriangle->C.y;
    else if( ptriangle->C.y < ymin) ymin = ptriangle->C.y;

    //set pointer
    pAabb=aabbs+i;

    //set values of aabb
    pAabb->min.x=xmin;
    pAabb->max.x=xmax;
    pAabb->min.y=ymin;
    pAabb->max.y=ymax;
    pAabb->min.z=ptriangle->A.z;
    pAabb->max.z=ptriangle->A.z+ptriangle->A.w;//w = height?
  }
}
