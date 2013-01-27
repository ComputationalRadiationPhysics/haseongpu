// Libraries
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "vector_types.h"
#include "assert.h"
#include "string.h"
#include <vector>
#include "curand_kernel.h"


// User header files
#include "datatypes.h"
#include "geometry.h"
#include "datatypes.h"
#include "generate_testdata.h"
#include "print.h"
#include "geometry_gpu.h"
#include "ase_bruteforce_kernel.h"
#include "ase_bruteforce_cpu.h"
#include "testdata_transposed.h"
#include "naive_ray_propagation.h"
#include "buildgrid.h"

int main(int argc, char **argv){
  const unsigned rays_per_sample = pow(2,18);
  const int threads = 256;
  char runmode[20];
  float runtime = 0.0;
  
  // Parse Commandline
  if(argc <= 1){
    fprintf(stderr, "C No commandline arguments found\n");
    fprintf(stderr, "C Usage    : ./octrace --mode=[runmode]\n");
    fprintf(stderr, "C Runmodes : bruteforce_gpu\n");
    fprintf(stderr, "             naive_ray_propagation\n");
    return 0;
  }
  
  // Generate testdata
  fprintf(stderr, "C Generate Testdata\n");
  std::vector<PrismCu>  *prisms = generatePrismsFromTestdata(host_mesh_z, host_p_in, host_size_p, host_t_in, host_size_t, host_mesh_z);
  std::vector<PointCu> *samples = generateSamplesFromTestdata(host_mesh_z, host_p_in, host_size_p);
  std::vector<double>    *betas = generateBetasFromTestdata(host_beta_v, host_mesh_z * host_size_t);
  std::vector<double>      *ase = new std::vector<double>(samples->size(), 0);
  const unsigned rays_total = rays_per_sample * samples->size();

  // Run 
  unsigned i;
  for(i=1; i < argc; ++i){
    if(strncmp(argv[i], "--mode=", 6) == 0){
      if(strstr(argv[i], "bruteforce_gpu") != 0){
  	runtime = runAseBruteforceGpu(samples, prisms, rays_per_sample, betas, ase, threads);
	strcpy(runmode, "Bruteforce GPU");

      }
      else if(strstr(argv[i], "naive_ray_propagation") != 0){
	runtime = runNaiveRayPropagation(ase);
	strcpy(runmode, "Naive Ray Propagation GPU");
	  }
	  else{
	fprintf(stderr, "C Runmode is not known\n");
	return 0;

      }

    }

  }

  // Print Solution
  unsigned sample_i;
  fprintf(stderr, "C Solutions\n");
  for(sample_i = 0; sample_i < ase->size(); ++sample_i){
    fprintf(stderr, "C ASE PHI of sample %d: %.80f\n", sample_i, ase->at(sample_i));

  }

  // Print statistics
  unsigned blocksPerSample = rays_per_sample / threads;
  unsigned blocks = blocksPerSample * samples->size();
  fprintf(stderr, "\n");
  fprintf(stderr, "C Statistics\n");
  fprintf(stderr, "C Prism             : %d\n", (int) prisms->size());
  fprintf(stderr, "C Triangles         : %d\n", (int) prisms->size() * 8);
  fprintf(stderr, "C Samples           : %d\n", (int) samples->size());
  fprintf(stderr, "C Rays/Sample       : %d\n", rays_per_sample);
  fprintf(stderr, "C GPU Blocks        : %d\n", blocks);
  fprintf(stderr, "C GPU Threads       : %d\n", threads);
  fprintf(stderr, "C GPU Blocks/Sample : %d\n", blocksPerSample);
  fprintf(stderr, "C Runmode           : %s \n", runmode);
  fprintf(stderr, "C Runtime           : %f s\n", runtime / 1000.0);
  fprintf(stderr, "\n");

  return 0;
}


