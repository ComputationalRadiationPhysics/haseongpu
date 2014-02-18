/**
 * @brief Propagates an ray through the triangle/prism/crystal structure.
 *        On each step the next triangle on the ray path will be calculated 
 *        from the current triangle (startTriangle in the beginning).
 *        length and startpoint of propagation is stored inside the
 *        ray struct. The propagation ends when the length of the ray 
 *        is reduced to zero. It is possible to do propagation with
 *        or without reflection. In case of reflection, the rays will
 *        have a longer way.
 * 
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @author Marius Melzer
 * @licence GPLv3
 *
 */

#ifndef propagate_ray_H
#define propagate_ray_H

#include <mesh.h>
#include <reflection.h> /* ReflectionPlane */

/**
 * @brief Direct ray propagation without reflection
 *
 * @param ray           The ray which will propagate through the prisms.
 * @param startLevel    The level where the startpoint of the ray is located.
 * @param startTriangle The triangle where the startpoint of the ray is locaed
 * @param mesh          All information about triangles, points, contants. 
 *                      See mesh.h for details.
 * @param sigmaA        Absorption value of the ray.
 * @param sigmaE        Emission value of the ray.
 *
 * @return gain         Integral of ray length multiplied by beta values.
 *                      See at the code or accompanying paper for more clarity.
 *
 */
__device__ double propagateRay(Ray ray, 
			       unsigned *startLevel, 
			       unsigned  *startTriangle, 
			       const Mesh &mesh,
			       const double sigmaA, 
			       const double sigmaE
			       );
/**
 * @brief Indirect ray propagation with reflections on upper and lower surface
 *
 * @param startPoint      Point where the ray should start from.
 * @param endPoint        Point where the will end.
 * @param reflections      Number of reflections the ray will do
 *                        from startPoint to endPoint.
 * @param reflectionPlane Plane of first reflection (upper or lower surface of gain medium)
 * @param startLevel      The level where the startpoint of the ray is located.
 * @param startTriangle   The triangle where the startpoint of the ray is locaed
 * @param mesh            All information about triangles, points, contants. 
 *                        See mesh.h for details.
 * @param sigmaA          Absorption value of the ray.
 * @param sigmaE          Emission value of the ray.
 *
 * @return gain           Integral of ray length multiplied by beta values.
 *                        See at the code or accompanying paper for more clarity.
 *
 */
__device__ double propagateRayWithReflection(Point startPoint, 
					     const Point endPoint, 
					     const unsigned reflections,
					     ReflectionPlane reflectionPlane,
					     unsigned startLevel, 
					     unsigned startTriangle, 
					     const Mesh &mesh, 
					     const double sigmaA, 
					     const double sigmaE);

#endif /* propagate_ray_H */


