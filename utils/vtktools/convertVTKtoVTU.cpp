#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <iostream>
#include <string>


int main(int argc,char* argv[]){
  if(argc < 2) {
    std::cerr << "needs a vtk argument!" << std::endl;
    return 1;
  }

  std::string fileName = argv[1];


  vtkUnstructuredGridReader *inputReader = vtkUnstructuredGridReader::New();
  inputReader->SetFileName(fileName.c_str());

  vtkXMLUnstructuredGridWriter *writer = vtkXMLUnstructuredGridWriter::New();

  char* pos = strtok(argv[1], ".");
  std::string basename = std::string(pos);

  writer->SetInput(inputReader->GetOutput());
  writer->SetFileName( (basename + ".vtu").c_str());
  writer->Write();

  return 0;

}
