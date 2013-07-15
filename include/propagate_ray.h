#ifndef propagate_ray_H
#define propagate_ray_H

#include <mesh.h>

__host__ __device__ double propagateRay(Ray ray, 
					unsigned startLevel, 
					Triangle startTriangle, 
					const double sigmaA, 
					const double sigmaE, 
					const double nTot, 
					const double thickness);

#endif /* propagate_ray_H */


