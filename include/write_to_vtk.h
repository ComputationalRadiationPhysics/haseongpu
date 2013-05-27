#include <iostream>
#include <fstream>


int writeToVtk(std::vector<float> * points, unsigned numberOfPoints){
  std::ofstream vtkFile;
  vtkFile.open("octrace.vtk");
  vtkFile << "# vtk DataFile Version 2.0" << std::endl;
  vtkFile << "octrace vtk file" << std::endl;
  vtkFile << "ASCII" << std::endl;
  vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
  vtkFile << "POINTS " << points->size() <<  "float" << std::endl;
  
  
  for(int i=0; i < points->size(); ++i){
    vtkFile << std::endl;

  }
  
  



}
