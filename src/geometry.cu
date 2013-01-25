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
				   PointCu rayObjective, //Richtungsvektor des Strahls
				   PointCu p1, //1.Punkt des Dreiecks
				   PointCu p2, //2.Punkt des Dreiecks
				   PointCu p3) //3.Punkt des Dreiecks
{
  double s2; //2.barizentrische Koordinate des Dreiecks
  double s3; //3.barizentrische Koordinate des Dreiecks
  //1.barizentrische Koordinate des Dreiecks ergibt sich mit 1.-s2-s3
  double t; //Geradenparameter
  //PointCu intersectionPoint = {0, 0, 0, 0};

  //Grenzwert fuer numerische Stabilitaet
  const double eps = 1e-6; //empirischer Wert, bei Moeller/Trumbore 1e-6

  //Variable fuer Determinante
  double determinante;

  //side12 und side13 sind Vektoren der Seiten
  //cross ist eine Hilfsvariable
  VectorCu side12, side13, rayDirection, cross;

  //Berechnung von Vektoren der Seiten:
  //1.Seite side12 von p1 nach p2
  side12.x = p2.x - p1.x;
  side12.y = p2.y - p1.y;
  side12.z = p2.z - p1.z;

  //2.Seite side13 von p1 nach p3
  side13.x = p3.x - p1.x;
  side13.y = p3.y - p1.y;
  side13.z = p3.z - p1.z;

  rayDirection.x = rayObjective.x - rayOrigin.x;
  rayDirection.y = rayObjective.y - rayOrigin.y;
  rayDirection.z = rayObjective.z - rayOrigin.z;

  //Gleichsetzen von Gereadengleichung und Ebenengleichung
  //Modell:  sp=p1+s2*(p2-p1)+s3*(p3-p1)=rayOrigin+t*rayDirection
  //Berechnung mit Cramerscher Regel zum Loesen des Gleichungssystems:
  //s2*(p2-p1)+s3*(p3-p1)-t*rayDirection = rayOrigin-p1 -> zu bestimmende Parameter: s2,s3,-t 

  //Kreuzprodukt von side13 und rayDirection
  cross.x = side13.y * rayDirection.z - side13.z * rayDirection.y;
  cross.y = side13.z * rayDirection.x - side13.x * rayDirection.z;
  cross.z = side13.x * rayDirection.y - side13.y * rayDirection.x;

  //Berechnung der Determinante mit Skalarprodukt
  determinante = cross.x * side12.x + cross.y * side12.y + cross.z * side12.z;


  //Test auf Parallelitaet
  //numerische Stabilitaet!!!

  if (determinante > -eps && determinante < eps){
    return -1.;
  }

  //Abstand Ursprung des Strahls zu p1
  VectorCu p1_rayOrigin; //=rayOrigin-p1;
  p1_rayOrigin.x = rayOrigin.x - p1.x;
  p1_rayOrigin.y = rayOrigin.y - p1.y;
  p1_rayOrigin.z = rayOrigin.z - p1.z;

  //barizentrische Koordinaten
  // sp=s1*p1+s2*p2+s3*p3
  //2. barizentrische Koordinate s2
  //=Skalarprodukt von p1_rayOrigin und cross
  s2 = cross.x * p1_rayOrigin.x + cross.y * p1_rayOrigin.y + cross.z * p1_rayOrigin.z;

  //Hilfsvariable
  VectorCu tempcross;
  //zunaenaehst Kreuzprodukt von rayDirection und side12
  tempcross.x = rayDirection.y * side12.z - rayDirection.z * side12.y;
  tempcross.y = rayDirection.z * side12.x - rayDirection.x * side12.z;
  tempcross.z = rayDirection.x * side12.y - rayDirection.y * side12.x;

  //s3=Skalarprodukt von rayDirection und side12
  //s3=(rayDirection x side12) *p1_rayOrigin
  s3 = tempcross.x * p1_rayOrigin.x + tempcross.y * p1_rayOrigin.y + tempcross.z * p1_rayOrigin.z;

  //Cramersche Regel -> Division durchfuehren
  double invdet = 1. / determinante;

  s2 = invdet*s2;
  s3 = invdet*s3;

  //weitere Verwendung der Hilfsvariable fuer Berechnung des Geradenparameters t
  //zunaechst Kreuzprodukt von side13 und side12
  tempcross.x = side13.y * side12.z - side13.z * side12.y;
  tempcross.y = side13.z * side12.x - side13.x * side12.z;
  tempcross.z = side13.x * side12.y - side13.y * side12.x;

  //t ist dann das Skalarprodukt von tempcross und p1_rayOrigin
  //t=(seite13,seite12) *p1_rayOrigin = -(seite12,seite13) *p1_rayOrigin
  t = tempcross.x * p1_rayOrigin.x + tempcross.y * p1_rayOrigin.y + tempcross.z * p1_rayOrigin.z;

  t = invdet*t;

  //Test,ob der Schnittpunkt innerhalb des Dreiecks liegt:

  //Ueberschereitungstest fuer barizentrische Koordinaten
  if (s2 < 0. || s2 > 1.) return -1.;

  //Ueberschereitungstest fuer barizentrische Koordinaten
  if (s3 < 0. || s3 > 1.) return -1.;

  //0 <= s1=1-s2-s3 <= 1 -> s2+s3<1   (s2+s3>0 schon durchgefuehrt,da s2>0 s3>0)
  if (s2 + s3 > 1.) return -1.;

  //Test, ob Strahl in Richtung des Dreiecks zeigt:
  // if (t < 0.) return 0;

  //Schnittpunktberechnung
  /* intersectionPoint.x = rayOrigin.x + t * rayDirection.x; */
  /* intersectionPoint.y = rayOrigin.y + t * rayDirection.y; */
  /* intersectionPoint.z = rayOrigin.z + t * rayDirection.z; */
  /* intersectionPoint.w = 1; */

  return t;
 
}

/**
   @attention slower than commented collide_prism_gpu but more pretty
**/
__device__ float collide_prism_gpu(PrismCu pr, RayCu r){
  //bool has_collide;
  /* PointCu first_intersection = {0, 0, 0, 0}; */
  /* PointCu intersection_point = {0, 0, 0, 0}; */
  double t1 = 0.;
  double t2 = 0.;
  PointCu A1 = pr.t1.A;
  PointCu B1 = pr.t1.B;
  PointCu C1 = pr.t1.C;
  PointCu A2 = {pr.t1.A.x, pr.t1.A.y, pr.t1.A.z + pr.t1.A.w, 1};
  PointCu B2 = {pr.t1.B.x, pr.t1.B.y, pr.t1.B.z + pr.t1.B.w, 1};
  PointCu C2 = {pr.t1.C.x, pr.t1.C.y, pr.t1.C.z + pr.t1.C.w, 1};
  float ray_distance = distance_gpu(r.P, r.direction);
  TriangleCu triangles[8] = {
    pr.t1,
    {A2, B2, C2},
    {A1, B1, A2},
    {B1, B2, A2},
    {B1, C1, C2},
    {B1, B2, C2},
    {A1, C1, C2},
    {A1, A2, C2}};

  float4 rayDirection;
  bool firstIntersectionFound = false;
  rayDirection.x = r.direction.x - r.P.x;
  rayDirection.y = r.direction.y - r.P.y;
  rayDirection.z = r.direction.z - r.P.z;
  float absRay = rayDirection.x * rayDirection.x + rayDirection.y * rayDirection.y + rayDirection.z * rayDirection.z;
  absRay = sqrtf(absRay);

  // test for collision on all triangles of an prism
  unsigned i; 
  for(i = 0; i < 8; ++i){
    t1 = intersectionRayTriangleGPU(r.P, r.direction, triangles[i].A, triangles[i].B, triangles[i].C);

    if(t1 >= 0.){
      if(!firstIntersectionFound){
	firstIntersectionFound= true;
	t2 = t1;
      }
      else{
	// Filter double collisions
	if(fabs(t2 - t1) > 0.0000001){
	   if(t1 > 1. && t2 > 1.)
	    return 0.;
	   if(t1>1.) t1=1.;
	   if(t2>1.) t2=1.;

	  return fabs(t2 - t1) * absRay; 
	}

      }

    }

  }

  return 0.;
}
