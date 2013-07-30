#include <iostream> /* cerr */
#include <fstream> /* ofstream */
#include <vector> /* vector */
#include <iomanip> /* std::setprecision() */
#include <mesh.h>
#include <cstdlib> /* atof */
#include <string>

int writeToVtk(Mesh *mesh,
	       std::vector<double>* ase,
	       std::string filename){

  std::cerr << "C Write experiment data to vtk-file" << std::endl;
  std::ofstream vtkFile;
  vtkFile.open(filename.c_str());

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

int compareVtk(std::vector<double> *ase, std::string filename){
  std::ifstream filestream;
  std::string line;
  bool foundLine = false;
  double value = 0;
  double diff = 0;
  unsigned i = 0;
  double minDiff = 10000; // should be enough
  double maxDiff = 0;
  double totalDiff = 0;
  double aseTotal = 0;
  double smallDiff = 0.1;
  

  // No compare vtk was given
  if(!filename.compare("")){
    return 0;
  }

  // Sum up ase values
  for(unsigned i = 0; i < ase->size(); ++i){
    aseTotal += ase->at(i);
  }

  filestream.open(filename.c_str(), std::ifstream::in);

  if(filestream.is_open()){
    while(filestream.good()){
      std::getline(filestream, line);
      std::size_t found = line.find("LOOKUP_TABLE default");
      if(found != std::string::npos){ 
	foundLine = true;
	std::getline(filestream, line);
      }
      if(foundLine){
	if(i == ase->size())
	  break;
	value = (double) atof(line.c_str());
	totalDiff += abs(ase->at(i) - value);
	diff = abs(ase->at(i) / value - 1);
	ase->at(i) = diff;

	if(diff >= maxDiff)
	  maxDiff = diff;

	if(diff <= minDiff)
	  minDiff = diff;

	if(diff >= smallDiff){
	  std::cerr << "C ASE relative difference[" << i << "]: " << diff << " > " << smallDiff << std::endl;
	}
	i++;

      }

    }

  }
  else{
    std::cerr << "C Can't open file " << filename << " for comparison" << std::endl;
    return 1;
  }

  std::cerr << "C ASE max. relative difference: " << maxDiff << std::endl;
  std::cerr << "C ASE min. relative difference: " << minDiff << std::endl;
  std::cerr << "C ASE tot. relative difference: " << totalDiff / aseTotal << std::endl;
  filestream.close();
  return 0;
}
