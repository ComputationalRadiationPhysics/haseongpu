#include <iostream> /* cerr */
#include <fstream> /* ofstream */
#include <vector> /* vector */


int writeToVtk(std::vector<int>* points, 
	       unsigned numberOfPoints,  
	       std::vector<unsigned>* triangleIndices, 
	       unsigned numberOfTriangles, 
	       unsigned numberOfLevels, 
	       std::vector<double>* ase){

  std::cerr << "C Write experiment data to vtk-file" << std::endl;
  std::ofstream vtkFile;
  vtkFile.open("octrace.vtk");

  // Write header of vtk file
  vtkFile << "# vtk DataFile Version 2.0" << std::endl;
  vtkFile << "octrace vtk file" << std::endl;
  vtkFile << "ASCII" << std::endl;

  // Write point data
  vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
  vtkFile << "POINTS " << points->size() <<  " float" << std::endl;
  for(unsigned level_i=0; level_i <= numberOfLevels; ++level_i){
    for(unsigned point_i=0; point_i < numberOfPoints; ++point_i){
      vtkFile << points->at(point_i) << " " << points->at(point_i + numberOfPoints) << " " << level_i << std::endl;

    }

  }

  // Write prism data
  vtkFile << "CELLS" << " " << (numberOfLevels-1) * numberOfTriangles << " " << (numberOfLevels-1) * numberOfTriangles * 7 << std::endl;
  for(unsigned level_i=0; level_i < numberOfLevels; ++level_i){
    for(unsigned triangle_i=0; triangle_i < numberOfTriangles; ++triangle_i){
      vtkFile << level_i * numberOfPoints + triangleIndices->at(triangle_i) << " "
	      << level_i * numberOfPoints + triangleIndices->at(numberOfTriangles + triangle_i) << " "
	      << level_i * numberOfPoints + triangleIndices->at(2 * numberOfTriangles + triangle_i) << " "
	      << level_i+1 * numberOfPoints + triangleIndices->at(triangle_i) << " "
	      << level_i+1 * numberOfPoints + triangleIndices->at(numberOfTriangles + triangle_i) << " "
	      << level_i+1 * numberOfPoints + triangleIndices->at(2 * numberOfTriangles + triangle_i) << std::endl;
	

    }

  }

  // Write cell type
  vtkFile << "CELL_TYPES " << (numberOfLevels-1) * numberOfTriangles << std::endl;
  for(unsigned prism_i=0; prism_i < (numberOfLevels-1) * numberOfTriangles; ++prism_i){
    vtkFile << "13" << std::endl;
  }

  // Write ase phi
  vtkFile << "POINT_DATA " << numberOfLevels * numberOfPoints << std::endl;
  vtkFile << "SCALARS scalars float 1 " << std::endl;
  vtkFile << "LOOKUP_TABLE default" << std::endl;

  for(unsigned ase_i=0; ase_i < numberOfLevels * numberOfPoints; ++ase_i){
    vtkFile << ase->at(ase_i) << std::endl;
  }
  
  vtkFile.close();

  return 0;
}
