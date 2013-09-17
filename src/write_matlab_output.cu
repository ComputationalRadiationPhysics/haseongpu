#include "write_matlab_output.h"
#include <fstream>
#include <vector>
#include <iomanip>

void writeMatlabOutput(
    const std::string experimentPath,
    const std::vector<float> ase,
    const std::vector<unsigned> N_rays, 
    const std::vector<double> expectedValues,
    const unsigned numberOfWavelengths,
    const unsigned numberOfSamples){

  std::ofstream aseFile;
  std::ofstream raysFile;
  std::ofstream expectedValuesFile;


  aseFile.open((experimentPath + "phi_ASE.txt").c_str());
  for(unsigned i = 0; i < numberOfSamples; ++i){
    for(unsigned j = 0; j < numberOfWavelengths; j++){
      aseFile << std::fixed << std::setprecision(20) << ase.at(i+j*numberOfSamples) << " ";
    }
    aseFile << std::endl;
  }
  aseFile.close();


  raysFile.open((experimentPath + "N_rays.txt").c_str());
  for(unsigned i = 0; i < numberOfSamples; ++i){
    for(unsigned j = 0; j<numberOfWavelengths; ++j){
      raysFile << N_rays.at(j) << " ";
      // TODO: change so that each samplepoint has its own number of rays! (adaptive number of rays!)
    }
    raysFile << std::endl;
  }
  raysFile.close();


  expectedValuesFile.open((experimentPath + "expected_values.txt").c_str());
  for(unsigned i = 0; i < numberOfSamples; ++i){
    for(unsigned j = 0; j < numberOfWavelengths; j++){
      expectedValuesFile << std::fixed << std::setprecision(20) <<  expectedValues.at(i+j*numberOfSamples) << " " ;
    }
    expectedValuesFile << std::endl;
  }
  expectedValuesFile.close();
}
