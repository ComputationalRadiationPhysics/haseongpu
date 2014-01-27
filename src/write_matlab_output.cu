#include "write_matlab_output.h"
#include <fstream>
#include <vector>
#include <iomanip>

template <typename T>
void write3dMatrix(
    const std::vector<T>& data, 
    std::ofstream &file,
    const unsigned rowCount,
    const unsigned columnCount,
    const unsigned pageCount
    ){
  
  unsigned elementsPerPage = rowCount*columnCount;
  file << rowCount << " " << columnCount << " " << pageCount << std::endl;
  for(unsigned page_i = 0; page_i < pageCount; ++page_i){
    for(unsigned j = 0; j < elementsPerPage ; ++j){
      file << std::fixed << std::setprecision(20) << data.at(j + page_i * elementsPerPage) << " , ";
    }
  }
  file << std::endl;

}


void writeMatlabOutput(
    const std::string experimentPath,
    const std::vector<float> ase,
    const std::vector<unsigned> N_rays, 
    const std::vector<double> expectedValues,
    const unsigned numberOfSamples,
    const unsigned numberOfLevels){

  std::ofstream aseFile;
  std::ofstream raysFile;
  std::ofstream expectedValuesFile;
  const unsigned samplesPerLevel = numberOfSamples/numberOfLevels;

  aseFile.open((experimentPath + "phi_ASE.txt").c_str());
  write3dMatrix(ase,aseFile,samplesPerLevel,numberOfLevels,1);
  aseFile.close();

  raysFile.open((experimentPath + "N_rays.txt").c_str());
  write3dMatrix(N_rays,raysFile,samplesPerLevel,numberOfLevels,1);
  raysFile.close();

  expectedValuesFile.open((experimentPath + "mse_values.txt").c_str());
  write3dMatrix(expectedValues,expectedValuesFile,samplesPerLevel,numberOfLevels,1);
  expectedValuesFile.close();
}
