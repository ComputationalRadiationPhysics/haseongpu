#include <fstream>
#include <vector>
#include <iomanip>

void writeMatlabOutput(
    std::vector<float>* ase,
    std::vector<unsigned> *N_rays, 
    std::vector<double> *expectedValues,
    unsigned numberOfWavelengths,
    unsigned numberOfSamples){
  std::ofstream aseFile;
  std::ofstream raysFile;
  std::ofstream expectedValuesFile;


  aseFile.open("phi_ASE.txt");
  for(unsigned i = 0; i < numberOfSamples; ++i){
    for(unsigned j = 0; j < numberOfWavelengths; j++){
      aseFile << std::fixed << std::setprecision(20) << ase->at(i+j*numberOfSamples) << " ";
    }
    aseFile << std::endl;
  }
  aseFile.close();


//  if(N_rays != NULL) {
	  raysFile.open("N_rays.txt");
	  for(unsigned i = 0; i < numberOfSamples; ++i){
		  for(unsigned j = 0; j<numberOfWavelengths; ++j){
			  raysFile << N_rays->at(j) << " ";
			  // TODO: change so that each samplepoint has its own number of rays! (adaptive number of rays!)
		  }
		  raysFile << std::endl;
	  }
	  raysFile.close();
 // }


  //if(expectedValues != NULL){
	  expectedValuesFile.open("expected_values.txt");
	  for(unsigned i = 0; i < numberOfSamples; ++i){
		  for(unsigned j = 0; j < numberOfWavelengths; j++){
			  expectedValuesFile << std::fixed << std::setprecision(20) <<  expectedValues->at(i+j*numberOfSamples) << " " ;
		  }
		  expectedValuesFile << std::endl;
	  }
	  expectedValuesFile.close();
 // }
}
