#include <iostream> /* cerr */
#include <fstream> /* ofstream */
#include <vector> /* vector */

#ifndef WRITE_TO_VTK_H
#define WRITE_TO_VTK_H

/**
 * @brief creates a VTK file based on the mesh structure and the dnd_ase values
 *
 * @param *points the sample points of the structure
 *
 * @param numberOfPoints the number of points in one level of the mesh
 *
 * @param *triangleIndices indices to the triangle coordinates
 *
 * @param numberOfTriangles the number of triangles in one layer of the mesh
 *
 * @param numberOfLevels the number of layers
 *
 * @param thicknessOfPrism the thickness of one layer
 *
 * @param *ase the dndt_ase values for each sample point
 *
 * @return 0
 *
 */

int writeToVtk(std::vector<double>* points, 
	       unsigned numberOfPoints,  
	       std::vector<unsigned>* triangleIndices, 
	       unsigned numberOfTriangles, 
	       unsigned numberOfLevels,
	       float thicknessOfPrism,
	       std::vector<double>* ase);

#endif /* WRITE_TO_VTK_H */
