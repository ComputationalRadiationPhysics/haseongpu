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
    const unsigned numberOfSamples,
    const unsigned numberOfLevels){

  std::ofstream aseFile;
  std::ofstream raysFile;
  std::ofstream expectedValuesFile;

  const unsigned samplesPerLevel = numberOfSamples/numberOfLevels;

  aseFile.open((experimentPath + "phi_ASE.txt").c_str());
  aseFile << samplesPerLevel << " " << numberOfLevels << " " << numberOfWavelengths << std::endl;
  for(unsigned wave_i = 0 ; wave_i < numberOfWavelengths; ++wave_i){
    for(unsigned ind = 0 ; ind < numberOfSamples; ++ind){
      aseFile << std::fixed << std::setprecision(20) << ase.at(ind + wave_i*numberOfSamples) << " , ";
    }
  }
  aseFile << std::endl;
  aseFile.close();


  raysFile.open((experimentPath + "N_rays.txt").c_str());
  raysFile << samplesPerLevel << " " << numberOfLevels << " " << numberOfWavelengths << std::endl;
  for(unsigned wave_i = 0 ; wave_i < numberOfWavelengths; ++wave_i){
    for(unsigned ind = 0 ; ind < numberOfSamples; ++ind){
      raysFile << N_rays.at(ind + wave_i*numberOfSamples) << " , ";
    }
    raysFile << std::endl;
  }
  raysFile.close();


  expectedValuesFile.open((experimentPath + "expected_values.txt").c_str());
  expectedValuesFile << samplesPerLevel << " " << numberOfLevels << " " << numberOfWavelengths << std::endl;
  for(unsigned wave_i = 0 ; wave_i < numberOfWavelengths; ++wave_i){
    for(unsigned ind = 0 ; ind < numberOfSamples; ++ind){
      expectedValuesFile << std::fixed << std::setprecision(20) <<  expectedValues.at(ind+wave_i*numberOfSamples) << " , " ;
    }
    expectedValuesFile << std::endl;
  }
  expectedValuesFile.close();



}
