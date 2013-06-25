#include "print_geometry.h"

void print_point(PointCu p){
  fprintf(stderr, "Point (%f, %f, %f, %f)\n",p.x, p.y, p.z, p.w);
}

void print_vector(VectorCu v){
  fprintf(stderr, "Vector (%f, %f, %f)\n",v.x, v.y, v.z);
}

void print_ray(RayCu r) {
  fprintf(stderr, "Ray (%f, %f, %f) -> (%f, %f, %f)\n",r.P.x, r.P.y, r.P.z, r.direction.x, r.direction.y, r.direction.z);
}

void print_plane(PlaneCu pl){
  fprintf(stderr, "Plane: \n\t");
  print_point(pl.P);
  fprintf(stderr, "\t");
  print_vector(pl.normal);
}

void print_triangle(TriangleCu tr){
  fprintf(stderr, "Triangle: \n");
  fprintf(stderr, "\t");  print_point(tr.A);
  fprintf(stderr, "\t");  print_point(tr.B);
  fprintf(stderr, "\t");  print_point(tr.C);
}
