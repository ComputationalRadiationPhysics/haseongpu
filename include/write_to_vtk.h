#include <iostream> /* cerr */
#include <fstream> /* ofstream */
#include <vector> /* vector */
#include <string>
#include <mesh.h>

#ifndef WRITE_TO_VTK_H
#define WRITE_TO_VTK_H

/**
 * @brief creates a VTK file based on the mesh structure and the dnd_ase values
 *
 * @param *mesh    all information about the mesh(points, triangles, constants etc.)
 *        *ase     the dndt_ase values for each sample point
 *        filename is the location where the vtk will be written
 *
 * @return 0
 *
 */
int writeToVtk(const Mesh& mesh,
	       const std::vector<double> ase,
	       const std::string filename, 
	       const unsigned raysPerSample,
	       const unsigned maxRaysPerSample,
	       const float expectationThreshold,
	       const bool useReflections,
	       const float runtime);

int writePrismToVtk(const Mesh& mesh,
	       const std::vector<double> prismData,
	       const std::string filename, 
	       const unsigned raysPerSample,
	       const unsigned maxRaysPerSample,
	       const float expectationThreshold,
	       const bool useReflections,
	       const float runtime);

std::vector<double> compareVtk(std::vector<double> ase, std::string filename, unsigned numberOfSamples);

#endif /* WRITE_TO_VTK_H */
