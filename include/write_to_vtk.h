#include <iostream> /* cerr */
#include <fstream> /* ofstream */
#include <vector> /* vector */

#ifndef WRITE_TO_VTK_H
#define WRITE_TO_VTK_H

int writeToVtk(std::vector<double>* points, 
	       unsigned numberOfPoints,  
	       std::vector<unsigned>* triangleIndices, 
	       unsigned numberOfTriangles, 
	       unsigned numberOfLevels,
	       float thicknessOfPrism,
	       std::vector<double>* ase);

#endif /* WRITE_TO_VTK_H */
