#include "stdio.h"

typedef struct point {
  float x;
  float y;
  float z;
} POINT;

typedef struct vector {
  float x;
  float y;
  float z;
} VECTOR;

typedef struct ray {
  point start;
  vector direction;
} RAY;

typedef struct triangle {
  point a;
  point b;
  point c;
} TRIANGLE;

typedef struct plane {
  point start;
  vector normal;
  
} PLANE;


point intersection(plane t, ray r);
void  print_point(point p);

int main(){

  triangle triangle_1 = {
    {0.0, 0.0, 0.0},
    {3.0, 0.0, 0.0},
    {0.0, 3.0, 0.0}};

  ray ray_1 = {
    {1.0, 1.0, 1.0},
    {0.0, 0.0, -1.0}};

  plane plane_1 = {
    {0.0, 0.0, 0.0},
    {0.0, 0.0, 1.0}};

  point point_1 = intersection(plane_1, ray_1);
  print_point(point_1);


  return 0;
}

/**
  @brief Intersection calculates the intersection between a plane p
         and a ray r. There is no detection for rays in the plane
	 or for parallel plane. 

  It uses the normal of the plane to derive the coordinate form 
  of the plane. With the help of a coordinate form it is very
  easy to get the intersection point between a ray and a plane.
**/
point intersection(plane p, ray r){
  point intersection_point = {0.0,0.0,0.0};
  float lampda = 0.0;
  float n1, n2, n3, x1, x2, x3, p1, p2, p3;

  // B has still to be calculated, not constant !
  float b = 1.0;

  n1 = p.normal.x;
  n2 = p.normal.y;
  n3 = p.normal.z;

  x1 = r.start.x;
  x2 = r.start.y;
  x3 = r.start.z;

  p1 = r.direction.x;
  p2 = r.direction.y;
  p3 = r.direction.z;

  lampda = (b - n1*x1 + n2*x2 * n3*x3) / (n1*p1 + n2*p2 + n3*p3);
  fprintf(stdout, "Lampda: %f\n", lampda);

  intersection_point.x = x1 + lampda * p1;
  intersection_point.y = x2 + lampda * p2;
  intersection_point.z = x3 + lampda * p3;

  return intersection_point;

}

void print_point(point p){
  fprintf(stdout, "Point\n");
  fprintf(stdout, "x: %f\n", p.x);
  fprintf(stdout, "y: %f\n", p.y);
  fprintf(stdout, "z: %f\n", p.z);

}
