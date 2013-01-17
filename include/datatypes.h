#ifndef DATATYPES_H
#define DATATYPES_H
typedef float4 point_cu;
typedef float4 vector_cu;

typedef struct triangle_cu{
  point_cu A;
  point_cu B;
  point_cu C;
} TRIANGLE_CU;

typedef struct prism_cu{
  triangle_cu t1; // height in w coordinate
} PRISM_CU;

typedef struct plane_cu {
  point_cu P;
  vector_cu normal;
} PLANE_CU;

typedef struct ray_cu {
  point_cu P;
  vector_cu direction;
} RAY_CU;

#endif /* DATATYPES_H */
