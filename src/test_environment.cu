#include <vector>
#include <mesh.h>
#include <thread>
#include <iostream>

void foo() 
{
  std::cout << "From Thread" << std::endl;
}


bool threadTest(){
  std::thread t0(foo);
  std::thread t1(foo);

}


bool testEnvironment (unsigned &threads, 
		       unsigned &blocks,
		       unsigned &hRaysPerSample,
		       const unsigned maxRaysPerSample,
		       const Mesh& dMesh,
		       const Mesh& hMesh,
		       const std::vector<double>& hSigmaA,
		       const std::vector<double>& hSigmaE,
		       const float expectationThreshold,
		       const bool useReflections,
		       std::vector<double> &dndtAse,
		       std::vector<float> &hPhiAse,
		       std::vector<double> &expectation
		       ){

  threadTest();

  return true;
}
