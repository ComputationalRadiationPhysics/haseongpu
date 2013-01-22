// Libraies
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "vector_types.h"
#include "assert.h"
#include "string.h"
#include <vector>

// User header files
#include "geometrie.h"
#include "datatypes.h"
#include "generate_testdata.h"
#include "print.h"
#include "geometry_gpu.h"
#include "ase_bruteforce_kernel.h"
#include "ase_bruteforce_cpu.h"

//----------------------------------------------------
// Host Code
//----------------------------------------------------
int main(int argc, char **argv){
  const unsigned max_rays = 256;
  const unsigned max_triangles = 2;
  const unsigned depth  = 1;
  const unsigned length = ceil(sqrt(max_triangles / 2));
  const int threads = 256;
  char* runmode = "";
  float runtime = 0.0;
  
  // Parse Commandline
  if(argc <= 1){
    fprintf(stderr, "C No commandline arguments found\n");
    fprintf(stderr, "C Usage    : ./octrace --mode=[runmode]\n");
    fprintf(stderr, "C Runmodes : bruteforce_cpu\n");
    fprintf(stderr, "             bruteforce_gpu\n");
    return 0;
  }
  
  // Generate testdata
  fprintf(stderr, "C Generate Testdata\n");
  std::vector<PrismCu> *prisms  = generate_prisms(length, length, depth);
  std::vector<PointCu> *samples = generate_samples(length, length, depth);
  std::vector<RayCu> *rays      = generate_sample_rays(length, length, depth, max_rays, samples);
  std::vector<float> *ase_cpu   = new std::vector<float>(samples->size(), 0);
  std::vector<float> *ase_gpu   = new std::vector<float>(samples->size(), 0);
  
  unsigned i;
  for(i=1; i < argc; ++i){
    if(strncmp(argv[i], "--mode=", 6) == 0){
      if(strstr(argv[i], "bruteforce_cpu") != 0){
  	runtime = run_ase_bruteforce_cpu(samples, prisms, rays, ase_cpu);
      }else
      if(strstr(argv[i], "bruteforce_gpu") != 0){
  	runtime = run_ase_bruteforce_gpu(samples, prisms, rays, ase_gpu, threads);
      }


    }
       
  }

  // Print statistics
  unsigned blocks_per_sample = ceil(rays->size() / threads);
  unsigned blocks = blocks_per_sample * samples->size();
  fprintf(stderr, "\n");
  fprintf(stderr, "C Prism             : %d\n", (int) prisms->size());
  fprintf(stderr, "C Triangles         : %d\n", (int) prisms->size() * 8);
  fprintf(stderr, "C Samples           : %d\n", (int) samples->size());
  fprintf(stderr, "C Rays/Sample       : %d\n", max_rays);
  fprintf(stderr, "C GPU Blocks        : %d\n", blocks);
  fprintf(stderr, "C GPU Threads       : %d\n", threads);
  fprintf(stderr, "C GPU Blocks/Sample : %d\n", blocks_per_sample);
  fprintf(stderr, "C Runmode           : %s \n", runmode);
  fprintf(stderr, "C Runtime           : %f s\n", runtime / 1000.0);
  fprintf(stderr, "\n");

  return 0;
}

