

#include "datatypes.h"


//p1 ... Punkt am rechten Winkel!!!!!
__device__ double intersectionRayRectGPU(PointCu rayOrigin, //Ursprung des Strahls
				   PointCu rayObjective, //Richtungsvektor des Strahls
				   PointCu p1, //1.Punkt des Rechtecks am rechten Winkel!!!
				   PointCu p2, //2.Punkt des Rechtecks
				   PointCu p3) //3.Punkt des Rechtecks
{
  double s2; //2.barizentrische Koordinate des Rechtecks
  double s3; //3.barizentrische Koordinate des Rechtecks

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


  //Test, ob Strahl in Richtung des Rechtecks zeigt:
  // if (t < 0.) return 0;

  //Schnittpunktberechnung
  /* intersectionPoint.x = rayOrigin.x + t * rayDirection.x; */
  /* intersectionPoint.y = rayOrigin.y + t * rayDirection.y; */
  /* intersectionPoint.z = rayOrigin.z + t * rayDirection.z; */
  /* intersectionPoint.w = 1; */

  return t;
 
}
