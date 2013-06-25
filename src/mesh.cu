#include "mesh.h"

float Ray::length() {
  return dir.length();
}

void Ray::normalize() {
  dir.normalize();
}

float Vector::length() {
  return sqrt(x*x + y*y + z*z);
}

void Vector::normalize() {
  int l = length();
  x = x/l;
  y = y/l;
  z = z/l;
}
