#ifndef propagate_ray_H
#define propagate_ray_H

#include <mesh.h>

/**
 * @brief Propagates an ray through the triangle/prism/crystal structure.
 *        On each step the next triangle on the ray path will be calculated 
 *        from the current triangle (startTriangle in the beginning).
 *        length and startpoint of propagation is stored inside the
 *        ray struct. The propagation ends when the length of the ray 
 *        is reduced to zero.
 * 
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @author Marius Melzer
 *
 * @licence GPLv3
 *
 * @param ray the ray which will propagate through the prisms
 *
 * @param startLevel the level where the startpoint of the ray is located
 *
 * @param startTriangle the triangle where the startpoint of the ray is locaed
 *
 * @param sigmaA
 * 
 * @param sigmaE
 *
 * @param nTot
 *
 * @param thickness is the thickness of one crystel slice (one level)
 *
 * @return gain Integral of ray length multiplied by beta values.
 *              See at the code or accompanying paper for more clarity.
 *
 **/
__host__ __device__ double propagateRay(Ray ray, 
					unsigned startLevel, 
					unsigned  startTriangle, 
					Mesh *mesh,
					const double sigmaA, 
					const double sigmaE
					);

#endif /* propagate_ray_H */


