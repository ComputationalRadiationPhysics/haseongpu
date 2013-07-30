#include <iostream> /* cerr */
#include <fstream> /* ofstream */
#include <vector> /* vector */
#include <iomanip> /* std::setprecision() */
#include <mesh.h>

int writeToVtk(Mesh *mesh,
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
  vtkFile << "POINTS " << mesh->numberOfSamples <<  " float" << std::endl;
  for(unsigned level_i=0; level_i < mesh->numberOfLevels; ++level_i){
    for(unsigned point_i=0; point_i < mesh->numberOfPoints; ++point_i){
      vtkFile << std::fixed << std::setprecision(6) << mesh->points[point_i] << " " << mesh->points[point_i + mesh->numberOfPoints] << " " << level_i * mesh->thickness << std::endl;

    }

  }

  // Write prism data
  vtkFile << "CELLS" << " " << mesh->numberOfPrisms << " " << mesh->numberOfPrisms * 7 << std::endl;
  for(unsigned level_i=0; level_i < (mesh->numberOfLevels - 1); ++level_i){
    for(unsigned triangle_i=0; triangle_i < mesh->numberOfTriangles; ++triangle_i){
      vtkFile << "6 " 
		  << level_i * mesh->numberOfPoints + mesh->triangles[triangle_i] << " "
	      << level_i * mesh->numberOfPoints + mesh->triangles[mesh->numberOfTriangles + triangle_i] << " "
	      << level_i * mesh->numberOfPoints + mesh->triangles[2 * mesh->numberOfTriangles + triangle_i] << " "
	      << (level_i+1) * mesh->numberOfPoints + mesh->triangles[triangle_i] << " "
	      << (level_i+1) * mesh->numberOfPoints + mesh->triangles[mesh->numberOfTriangles + triangle_i] << " "
	      << (level_i+1) * mesh->numberOfPoints + mesh->triangles[2 * mesh->numberOfTriangles + triangle_i] << std::endl;
	
    }

  }

  // Write cell type
  vtkFile << "CELL_TYPES " << mesh->numberOfPrisms << std::endl;
  for(unsigned prism_i=0; prism_i < mesh->numberOfPrisms; ++prism_i){
    vtkFile << "13" << std::endl;
  }

  // Write ase phi
  vtkFile << "POINT_DATA " << mesh->numberOfSamples << std::endl;
  vtkFile << "SCALARS scalars float 1" << std::endl;
  vtkFile << "LOOKUP_TABLE default" << std::endl;

  for(unsigned ase_i=0; ase_i < mesh->numberOfSamples; ++ase_i){
    vtkFile << std::fixed << std::setprecision(6) << ase->at(ase_i) << std::endl;
  }
  
  vtkFile.close();

  return 0;
}
