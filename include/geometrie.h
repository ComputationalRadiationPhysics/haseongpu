#include "datatypes.h"

#ifndef GEOMETRIE_H
#define GEOMETRIE_H
PointCu  collide_triangle(TriangleCu t, RayCu r);
float  collide_prism(PrismCu pr, RayCu r);
float4 to_barycentric(TriangleCu t, RayCu r);
PointCu intersection(PlaneCu p, RayCu r);
float distance(PointCu a, PointCu b);
VectorCu crossproduct(VectorCu a, VectorCu b);
float skalar_mul(VectorCu a, VectorCu b);

/**   
   @brief Calculates the barycentric coordinates of the triangle
          and the intersectionpoint of the ray and the triangle.
	  Algorithm based on a paper. exact copy of pseudo code.

   @return PointCu {0,0,0,0} if there is no intersection triangle/ray
   @return PointCu {x,y,z,1} barycentric coordinates of intersection triangle/ray
 **/
float4 to_barycentric(TriangleCu tr, RayCu ray){
  float4 b = {0,0,0,0};
  VectorCu e1, e2, q, s, r, ray_direction;
  PointCu p0, p1, p2;
  float a, f, u, v, t;


  p0 = tr.A;
  p1 = tr.B;
  p2 = tr.C;
  
  ray_direction.x = ray.direction.x - ray.P.x;
  ray_direction.y = ray.direction.y - ray.P.y;
  ray_direction.z = ray.direction.z - ray.P.z;

  e1.x = p1.x - p0.x;
  e1.y = p1.y - p0.y;
  e1.z = p1.z - p0.z;

  e2.x = p2.x - p0.x;
  e2.y = p2.y - p0.y;
  e2.z = p2.z - p0.z;

  q = crossproduct(ray_direction, e2);
  a = skalar_mul(e1, q);
  
  // a is to close to 0
  if(fabs(a) < 0.000001)
    return b;

  f = 1 / a;
  
  s.x = ray.P.x - p0.x;
  s.y = ray.P.y - p0.y;
  s.z = ray.P.z - p0.z;

  u = f * skalar_mul(s, q);

  if(u < 0.0)
    return b;

  r = crossproduct(s, e1);
  v = f * skalar_mul(ray_direction, r);
  if( v < 0.0 || (u + v) > 1)
    return b;
  
  t = f * skalar_mul(e2, q);
  
  b.x = u;
  b.y = v;
  b.z = t;
  b.w = 1;

  return b;
}

/**
   @brief Detects collisions of a triangle and a ray without
   a precondition.
   
   @return PointCu {0,0,0,0} if there is no intersection triangle/ray
   @return PointCu {x,y,z,1} barycentric coordinates of intersection triangle/ray
**/
PointCu collide_triangle(TriangleCu t, RayCu r){
  PlaneCu pl;
  float b1, b2, b3, c1, c2, c3;

  b1 = t.B.x - t.A.x;
  b2 = t.B.y - t.A.y;
  b3 = t.B.z - t.A.z;

  c1 = t.C.x - t.A.x;
  c2 = t.C.y - t.A.y;
  c3 = t.C.z - t.A.z;

  pl.P = t.A;
  pl.normal.x = (b2*c3 - b3*c2);
  pl.normal.y = (b3*c1 - b1*c3);
  pl.normal.z = (b1*c2 - b2*c1);

  float4 b = to_barycentric(t, r);
  // Maybe we can calculate intersection be barycentric coords
  PointCu p = intersection(pl, r);
  if(b.w == 1 && p.w == 1){
    return p;
  }
  else{
    PointCu no_inter = {0,0,0,0}; 
    return no_inter;
  }

}
/**
   @brief Detects collisions of a prism and ray
   
   @return float intersection distance
   @return 0 if in case of no intersection
**/
float collide_prism(PrismCu pr, RayCu r){
  //bool has_collide;
  PointCu first_intersection = {0, 0, 0, 0};
  PointCu intersection_point = {0, 0, 0, 0};
  PointCu A1 = pr.t1.A;
  PointCu B1 = pr.t1.B;
  PointCu C1 = pr.t1.C;
  PointCu A2 = {pr.t1.A.x, pr.t1.A.y, pr.t1.A.z + pr.t1.A.w, 1};
  PointCu B2 = {pr.t1.B.x, pr.t1.B.y, pr.t1.B.z + pr.t1.B.w, 1};
  PointCu C2 = {pr.t1.C.x, pr.t1.C.y, pr.t1.C.z + pr.t1.C.w, 1};
  float ray_distance = distance(r.P, r.direction);

  TriangleCu triangles[8] = {
    pr.t1,
    {A2, B2, C2},
    {A1, B1, A2},
    {B1, B2, A2},
    {B1, C1, C2},
    {B1, B2, C2},
    {A1, C1, C2},
    {A1, A2, C2}};


  // test for collision on all triangles of an prism
  unsigned i; 
  for(i = 0; i < 8; ++i){
    intersection_point = collide_triangle(triangles[i], r);
    if(intersection_point.w != 0){
      if(first_intersection.w == 0){
	first_intersection = intersection_point;
      }
      else{
	// Filter double Collision on edges or vertices
	if(first_intersection.x != intersection_point.x || first_intersection.y != intersection_point.y || first_intersection.z != intersection_point.z){

	  /*
	  if(distance(r.P, first_intersection) <= ray_distance && distance(r.P, intersection_point) > ray_distance)
	    return distance(r.direction, first_intersection);

	  if(distance(r.P, first_intersection) >= ray_distance && distance(r.P, intersection_point) < ray_distance)
	    return distance(r.direction, intersection_point);
	  
	  if(distance(r.P, first_intersection) > ray_distance || distance(r.direction, first_intersection) > ray_distance)
	    return 0;
	  */
	  return distance(first_intersection, intersection_point);
	}
      }

    }

  }

  return 0;
}

/**
   @brief Intersection calculates the intersection between a plane p
   and a ray r.

   It uses the normal of the plane to derive the coordinate form 
   of the plane. With the help of a coordinate form it is very
   easy to get the intersection point between a ray and a plane.

   ray   g: y~ = x~ + t*p~
   plane E: y~ = a~ + r*b~ + s*c~
   d  = n1*(x1+t*p1) + n2*(x2+t*p2) + n3*(x3+t*p3)
   d  = n~ * a~
**/
PointCu intersection(PlaneCu pl, RayCu r){
  PointCu intersection_point = {0, 0, 0, 0};

  float t, d;

  // vector coordinates
  float n1, n2, n3, x1, x2, x3, p1, p2, p3, a1, a2, a3;
  
  // just get the coordinates from the structs
  n1 = pl.normal.x;
  n2 = pl.normal.y;
  n3 = pl.normal.z;

  a1 = pl.P.x;
  a2 = pl.P.y;
  a3 = pl.P.z;

  x1 = r.P.x;
  x2 = r.P.y;
  x3 = r.P.z;

  p1 = r.direction.x - r.P.x;
  p2 = r.direction.y - r.P.y;
  p3 = r.direction.z - r.P.z;

  // calculation of intersection
  // this case for parallel rays, will be ignored for easier calculations
  float denominator = (n1*p1 + n2*p2 + n3*p3);
  if(abs(denominator) <= 0.000001)
    return intersection_point;

  d = n1*a1 + n2*a2 + n3*a3;
  t = (d - n1*x1 - n2*x2 - n3*x3) / denominator;

  // ignore intersections before the ray 
  if(t < 0)
  return intersection_point;

  intersection_point.x = x1 + t * p1;
  intersection_point.y = x2 + t * p2;
  intersection_point.z = x3 + t * p3;
  intersection_point.w = 1;
  return intersection_point;

}

float distance(PointCu a, PointCu b){
  float d = sqrt(pow((b.x - a.x), 2) + pow((b.y - a.y),2) + pow((b.z - a.z),2));
  return fabs(d);
}

VectorCu crossproduct(VectorCu a, VectorCu b){
  VectorCu c = {
    a.y*b.z - a.z*b.y,
    a.z*b.x - a.x*b.z,
    a.x*b.y - a.y*b.x
  };
  return c;
}

float skalar_mul(VectorCu a, VectorCu b){
  return a.x*b.x + a.y*b.y + a.z*b.z;
}
#endif /* GEOMETRIE_H */
