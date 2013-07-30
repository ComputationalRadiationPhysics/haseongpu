#include <iostream> /* cerr */
#include <fstream> /* ofstream */
#include <vector> /* vector */
#include <mesh.h>

#ifndef WRITE_TO_VTK_H
#define WRITE_TO_VTK_H

/**
 * @brief creates a VTK file based on the mesh structure and the dnd_ase values
 *
 * @param *mesh all information about the mesh(points, triangles, constants etc.)
 *        *ase  the dndt_ase values for each sample point
 *
 * @return 0
 *
 */

int writeToVtk(Mesh *mesh,
	       std::vector<double>* ase);

#endif /* WRITE_TO_VTK_H */
