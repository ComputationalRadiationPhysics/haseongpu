#include <fstream>
#include <vector>
#include <iomanip>

void writeDndtAse(std::vector<float>* ase){
  std::ofstream aseFile;
  aseFile.open("dndt_ASE.txt");

  for(unsigned i = 0; i < ase->size(); ++i){
    aseFile << std::fixed << std::setprecision(20) << ase->at(i) << std::endl;
  }

  aseFile.close();
}
