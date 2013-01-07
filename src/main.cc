#include "stdio.h"
#include "math.h"

//----------------------------------------------------
// Structures
//----------------------------------------------------
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

//----------------------------------------------------
// Auxillary function declaration
//----------------------------------------------------
point intersection(plane p, ray r);
bool  collide(triangle t, ray r);
bool  collide(triangle t, point p);
float distance(point a, point b);
void  print_point(point p);

//----------------------------------------------------
//  Calculations
//----------------------------------------------------
int main(){

  triangle triangle_1 = {
    {0.0, 0.0, 0.0},
    {3.0, 0.0, 0.0},
    {0.0, 3.0, 0.0}};

  ray ray_1 = {
    {5.0, 5.0, 5.0},
    {0.0, 0.0, 1.0}};

  plane plane_1 = {
    {0.0, 0.0, 0.0},
    {0.0, 0.0, 1.0}};

  if(collide(triangle_1, ray_1))
    fprintf(stdout, "Collision\n");
  else
    fprintf(stdout, "no collision\n");

  return 0;
}

//----------------------------------------------------
// Auxillary function definition
//----------------------------------------------------
/**
   @brief Detects collisions of triangle and point with
          precondition, that the point is on the same 
	  plane as the point.
 **/
bool collide(triangle t, point p){
  bool has_collide = true;
  float max_d = distance(t.a, t.b);
    if(distance(p, t.a) > max_d)
    has_collide = false;
  if(distance(p, t.b) > max_d) 
    has_collide = false;
  if(distance(p, t.c) > max_d) 
    has_collide = false;

  return has_collide;
}

/**
  @brief Detects collisions of a triangle and a ray without
         a precondition.
 **/
bool collide(triangle t, ray r){
  bool has_collide = true;
  plane pl;
  float b1, b2, b3, c1, c2, c3;

  b1 = t.b.x;
  b2 = t.b.y;
  b3 = t.b.z;

  c1 = t.c.x;
  c2 = t.c.y;
  c3 = t.c.z;

  pl.start = t.a;
  pl.normal.x = (b2*c3 - b3*c2);
  pl.normal.y = (b3*c1 - b1*c3);
  pl.normal.z = (b1*c2 - b2*c1);

  return collide(t, intersection(pl, r));
}

/**
  @brief Intersection calculates the intersection between a plane p
         and a ray r. There is no detection for rays in the plane
	 or for parallel plane. 

  It uses the normal of the plane to derive the coordinate form 
  of the plane. With the help of a coordinate form it is very
  easy to get the intersection point between a ray and a plane.

  ray   g: y~ = x~ + t*p~
  plane E: y~ = a~ + r*b~ + s*c~
           d  = n1*(x1+t*p1) + n2*(x2+t*p2) + n3*(x3+t*p3)
           d  = n~ * a~
**/
point intersection(plane pl, ray r){
  point intersection_point = {0.0,0.0,0.0};

  float t, d;

  // vector coordinates
  float n1, n2, n3, x1, x2, x3, p1, p2, p3, a1, a2, a3;
  
  // just get the coordinates from the structs
  n1 = pl.normal.x;
  n2 = pl.normal.y;
  n3 = pl.normal.z;

  a1 = pl.start.x;
  a2 = pl.start.y;
  a3 = pl.start.z;

  x1 = r.start.x;
  x2 = r.start.y;
  x3 = r.start.z;

  p1 = r.direction.x;
  p2 = r.direction.y;
  p3 = r.direction.z;

  // calculation of intersection
  d = n1*a1 + n2*a2 + n3*a3;
  t = (d - n1*x1 - n2*x2 - n3*x3) / (n1*p1 + n2*p2 + n3*p3);

  intersection_point.x = x1 + t * p1;
  intersection_point.y = x2 + t * p2;
  intersection_point.z = x3 + t * p3;

  return intersection_point;

}

float distance(point a, point b){
  float d = sqrt(pow((b.x - a.x), 2) + pow((b.y - a.y),2) + pow((b.z - a.z),2));
  return fabs(d);
}

void print_point(point p){
  fprintf(stdout, "Point\n");
  fprintf(stdout, "x: %f\n", p.x);
  fprintf(stdout, "y: %f\n", p.y);
  fprintf(stdout, "z: %f\n", p.z);

}
