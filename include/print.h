#ifndef PRINT_H
#define PRINT_H

void print_point(PointCu p);
void print_vector(VectorCu v);
void print_plane(PlaneCu pl);
void print_triangle(TriangleCu tr);


void print_point(PointCu p){
  fprintf(stderr, "Point (%f, %f, %f)\n",p.x, p.y, p.z);

}

void print_vector(VectorCu v){
  fprintf(stderr, "Vector (%f, %f, %f)\n",v.x, v.y, v.z);

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

#endif /* PRINT_H */
