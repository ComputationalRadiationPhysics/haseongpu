#include "datatypes.h"
#include "geometry_gpu.h"
#include "curand_kernel.h"

__device__ PointCu addVectorToPoint(PointCu p, VectorCu v) {
  PointCu result;

  result.x = p.x + v.x;
  result.y = p.y + v.y;
  result.z = p.z + v.z;

  return result;
}

__device__ RayCu generateRayGpu(PointCu sample, PrismCu startPrism, curandState localState){
  float u = curand_uniform(&localState);
  float v = curand_uniform(&localState);
  if((u+v) > 1){ //OPTIMIZE: remove if
    u = 1-u;
    v = 1-v;
  }
  const float w = 1-(u+v);

  PointCu A = startPrism.t1.A;
  PointCu B = startPrism.t1.B;
  PointCu C = startPrism.t1.C;

  // Get x and y coordinates from the random barycentric values
  const float xRand = u*A.x + v*B.x + w*C.x ;
  const float yRand = u*A.y + v*B.y + w*C.y ;

  // Take one of the given z-coordinates and add a random part of the prism height
  const float zRand = A.z + curand_uniform(&localState) * startPrism.t1.A.w;

  float ase=0.f;

  // Take the values to assemble a ray
  RayCu r = {
    sample,
    {xRand, yRand, zRand, ase}};
  return r;
}    

__device__ PrismCu selectPrism(int gid, PrismCu *prisms, int totalNumberOfPrisms){
  int totalNumberOfThreads = blockDim.x * gridDim.x;
  int threadsPerPrism = ceil( float(totalNumberOfThreads) / float(totalNumberOfPrisms) );
  int prism = gid / threadsPerPrism;

  return prisms[prism];
}   

__device__ float distance_gpu(PointCu a, PointCu b){
  float d = sqrt(pow((b.x - a.x), 2) + pow((b.y - a.y),2) + pow((b.z - a.z),2));
  return fabs(d);
}
__device__ VectorCu crossproduct_gpu(VectorCu a, VectorCu b){
  VectorCu c = {
    a.y*b.z - a.z*b.y,
    a.z*b.x - a.x*b.z,
    a.x*b.y - a.y*b.x
  };
  return c;
}


__device__ float skalar_mul_gpu(VectorCu a, VectorCu b){
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

__device__ double intersectionRayTriangleGPU(PointCu rayOrigin, //Ursprung des Strahls
    VectorCu rayDirection,
    PointCu p1, //1.Punkt des Dreiecks
    PointCu p2, //2.Punkt des Dreiecks
    PointCu p3) //3.Punkt des Dreiecks
{
  double s2; //2.barizentrische Koordinate des Dreiecks
  double s3; //3.barizentrische Koordinate des Dreiecks
  //1.barizentrische Koordinate des Dreiecks ergibt sich mit 1.-s2-s3
  double t = 0.; //Geradenparameter
  //PointCu intersectionPoint = {0, 0, 0, 0};

  //Grenzwert fuer numerische Stabilitaet
  const double eps = 1e-6; //empirischer Wert, bei Moeller/Trumbore 1e-6

  //Variable fuer Determinante
  double determinante;

  //side12 und side13 sind Vektoren der Seiten
  //cross ist eine Hilfsvariable
  VectorCu side12, side13, cross, p1_rayOrigin, tempcross;

  //Berechnung von Vektoren der Seiten:
  //1.Seite side12 von p1 nach p2
  side12.x = p2.x - p1.x;
  side12.y = p2.y - p1.y;
  side12.z = p2.z - p1.z;

  //2.Seite side13 von p1 nach p3
  side13.x = p3.x - p1.x;
  side13.y = p3.y - p1.y;
  side13.z = p3.z - p1.z;

  //Gleichsetzen von Gereadengleichung und Ebenengleichung
  //Modell:  sp=p1+s2*(p2-p1)+s3*(p3-p1)=rayOrigin+t*rayDirection
  //Berechnung mit Cramerscher Regel zum Loesen des Gleichungssystems:
  //s2*(p2-p1)+s3*(p3-p1)-t*rayDirection = rayOrigin-p1 -> zu bestimmende Parameter: s2,s3,-t 

  //Kreuzprodukt von side13 und rayDirection
  cross = crossproduct_gpu(side13, rayDirection);

  //Berechnung der Determinante mit Skalarprodukt
  determinante = skalar_mul_gpu(cross, side12);

  //Test auf Parallelitaet
  //numerische Stabilitaet!!!

  if (determinante > -eps && determinante < eps){
    return -1.;
  }

  //Abstand Ursprung des Strahls zu p1

  p1_rayOrigin.x = rayOrigin.x - p1.x;
  p1_rayOrigin.y = rayOrigin.y - p1.y;
  p1_rayOrigin.z = rayOrigin.z - p1.z;

  //barizentrische Koordinaten
  // sp=s1*p1+s2*p2+s3*p3
  //2. barizentrische Koordinate s2
  //=Skalarprodukt von p1_rayOrigin und cross
  s2 = skalar_mul_gpu(cross, p1_rayOrigin);

  //zunaenaehst Kreuzprodukt von rayDirection und side12
  tempcross = crossproduct_gpu(rayDirection, side12);

  //s3=Skalarprodukt von rayDirection und side12
  //s3=(rayDirection x side12) *p1_rayOrigin
  s3 = skalar_mul_gpu(tempcross, p1_rayOrigin);

  //Cramersche Regel -> Division durchfuehren
  double invdet = 1. / determinante;

  s2 = invdet*s2;
  s3 = invdet*s3;

  //Test,ob der Schnittpunkt innerhalb des Dreiecks liegt:
  //Ueberschereitungstest fuer barizentrische Koordinaten
  //if (s2 < 0. || s2 > 1. || s3 < 0. || s3 > 1. || s2 + s3 > 1.) return -1.;
  //s2 > 1. und s3 > 1. sind bereits enthalten in s2 + s3 > 1.
  if (s2 < 0. || s3 < 0. || s2 + s3 > 1.) return -1.;

  //weitere Verwendung der Hilfsvariable fuer Berechnung des Geradenparameters t
  //zunaechst Kreuzprodukt von side13 und side12
  tempcross = crossproduct_gpu(side13, side12);

  //t ist dann das Skalarprodukt von tempcross und p1_rayOrigin
  //t=(seite13,seite12) *p1_rayOrigin = -(seite12,seite13) *p1_rayOrigin
  t = skalar_mul_gpu(tempcross, p1_rayOrigin);
  t = invdet*t;

  return t;

}

__device__ double intersectionRayRectangleGPU(PointCu rayOrigin, //Ursprung des Strahls
    VectorCu rayDirection,
    PointCu p1, //1.Punkt des Dreiecks
    PointCu p2, //2.Punkt des Dreiecks
    PointCu p3) //3.Punkt des Dreiecks
{
  double s2; //2.barizentrische Koordinate des Dreiecks
  double s3; //3.barizentrische Koordinate des Dreiecks
  //1.barizentrische Koordinate des Dreiecks ergibt sich mit 1.-s2-s3
  double t = 0.; //Geradenparameter
  //PointCu intersectionPoint = {0, 0, 0, 0};

  //Grenzwert fuer numerische Stabilitaet
  const double eps = 1e-6; //empirischer Wert, bei Moeller/Trumbore 1e-6

  //Variable fuer Determinante
  double determinante;

  //side12 und side13 sind Vektoren der Seiten
  //cross ist eine Hilfsvariable
  VectorCu side12, side13, cross, p1_rayOrigin, tempcross;

  //Berechnung von Vektoren der Seiten:
  //1.Seite side12 von p1 nach p2
  side12.x = p2.x - p1.x;
  side12.y = p2.y - p1.y;
  side12.z = p2.z - p1.z;

  //2.Seite side13 von p1 nach p3
  side13.x = p3.x - p1.x;
  side13.y = p3.y - p1.y;
  side13.z = p3.z - p1.z;

  //Gleichsetzen von Gereadengleichung und Ebenengleichung
  //Modell:  sp=p1+s2*(p2-p1)+s3*(p3-p1)=rayOrigin+t*rayDirection
  //Berechnung mit Cramerscher Regel zum Loesen des Gleichungssystems:
  //s2*(p2-p1)+s3*(p3-p1)-t*rayDirection = rayOrigin-p1 -> zu bestimmende Parameter: s2,s3,-t 

  //Kreuzprodukt von side13 und rayDirection
  cross = crossproduct_gpu(side13, rayDirection);

  //Berechnung der Determinante mit Skalarprodukt
  determinante = skalar_mul_gpu(cross, side12);

  //Test auf Parallelitaet
  //numerische Stabilitaet!!!

  if (determinante > -eps && determinante < eps){
    return -1.;
  }

  //Abstand Ursprung des Strahls zu p1

  p1_rayOrigin.x = rayOrigin.x - p1.x;
  p1_rayOrigin.y = rayOrigin.y - p1.y;
  p1_rayOrigin.z = rayOrigin.z - p1.z;

  //barizentrische Koordinaten
  // sp=s1*p1+s2*p2+s3*p3
  //2. barizentrische Koordinate s2
  //=Skalarprodukt von p1_rayOrigin und cross
  s2 = skalar_mul_gpu(cross, p1_rayOrigin);

  //zunaenaehst Kreuzprodukt von rayDirection und side12
  tempcross = crossproduct_gpu(rayDirection, side12);

  //s3=Skalarprodukt von rayDirection und side12
  //s3=(rayDirection x side12) *p1_rayOrigin
  s3 = skalar_mul_gpu(tempcross, p1_rayOrigin);

  //Cramersche Regel -> Division durchfuehren
  double invdet = 1. / determinante;

  s2 = invdet*s2;
  s3 = invdet*s3;

  //Test,ob der Schnittpunkt innerhalb des Dreiecks liegt:
  if (s2 < 0. || s2 > 1. || s3 < 0. || s3 > 1.) return -1.;

  //weitere Verwendung der Hilfsvariable fuer Berechnung des Geradenparameters t
  //zunaechst Kreuzprodukt von side13 und side12
  tempcross = crossproduct_gpu(side13, side12);

  //t ist dann das Skalarprodukt von tempcross und p1_rayOrigin
  //t=(seite13,seite12) *p1_rayOrigin = -(seite12,seite13) *p1_rayOrigin
  t = skalar_mul_gpu(tempcross, p1_rayOrigin);
  t = invdet*t;

  return t;

}


__device__ float collide_prism_gpu(PrismCu pr, RayCu r, VectorCu rayDirection, double absRayDistance){
  PointCu A1 = pr.t1.A;
  PointCu B1 = pr.t1.B;
  PointCu C1 = pr.t1.C;
  PointCu A2 = {pr.t1.A.x, pr.t1.A.y, pr.t1.A.z + pr.t1.A.w, 1};
  PointCu B2 = {pr.t1.B.x, pr.t1.B.y, pr.t1.B.z + pr.t1.B.w, 1};
  PointCu C2 = {pr.t1.C.x, pr.t1.C.y, pr.t1.C.z + pr.t1.C.w, 1};
  double t[5];

  t[0] = intersectionRayTriangleGPU(r.P, rayDirection, A1, B1, C1);
  t[1] = intersectionRayTriangleGPU(r.P, rayDirection, A2, B2, C2);
  t[2] = intersectionRayRectangleGPU(r.P, rayDirection, A1, C1, A2);
  t[3] = intersectionRayRectangleGPU(r.P, rayDirection, A1, B1, A2);
  t[4] = intersectionRayRectangleGPU(r.P, rayDirection, C1, B1, C2);

  bool firstIntersectionFound = false;
  // test for collision on all triangles of an prism
  unsigned i; 
  for(i = 0; i < 5; ++i){
    if(t[i] >= 0.){
      if(!firstIntersectionFound){
        firstIntersectionFound= true;
        t[0] = t[i];
      }
      else{
        // Filter double collisions
        if(fabs(t[0] - t[i]) > 0.0000001){
          if(t[i] > 1. && t[0] > 1.)
            return 0.;
          if(t[i]>1.) t[i] = 1.;
          if(t[0]>1.) t[0] = 1.;

          return fabs(t[0] - t[i]) * absRayDistance; 
        }

      }

    }

  }

  return 0.;
}
