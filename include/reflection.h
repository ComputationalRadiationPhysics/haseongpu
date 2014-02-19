/**
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @licence GPLv3
 *
 */

#ifndef Reflection_H
#define Reflection_H

#include <mesh.h>
#include <geometry.h>

/**
 * @brief Calculates the reflectionPoint and reflectionAngle with upper or lower surface
 *        of gain medium. Depending on the number of reflection the
 *        ray still has to do and the reflection surface. 
 *
 * @param startPoint      Point where the ray should start from.
 * @param endPoint        Point where the will end.
 * @param reflectionsLeft Number of reflections the ray will do
 *                        from startPoint to endPoint.
 * @param reflectionPlane Plane of first reflection (upper or lower surface of gain medium).
 * @param reflectionPoint Pointer to the point the ray will be reflected.
 * @param reflectionAngle Pointer to the value of the angle the ray will be reflected with.
 *
 * @return 0 if everything goes right
 *         1 otherwise
 *
 */
__device__ int calcNextReflection(const Point startPoint, 
				  const Point endPoint, 
				  const unsigned reflectionsLeft, 
				  const ReflectionPlane reflectionPlane, 
				  Point *reflectionPoint, 
				  double *reflectionAngle, 
				  const Mesh &mesh);

#endif /* Reflection_H */
