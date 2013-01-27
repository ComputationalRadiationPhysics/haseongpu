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
  unsigned rays_total;
  char runmode[100];
  float runtime = 0.0;
  unsigned blocks = 0;
  unsigned threads = 0;
  // Parse Commandline
  if(argc <= 1){
    fprintf(stderr, "C No commandline arguments found\n");
    fprintf(stderr, "C Usage    : ./octrace --mode=[runmode] --rays=[number of rays]\n");
    fprintf(stderr, "C Runmodes : bruteforce_gpu\n");
    fprintf(stderr, "             naive_ray_propagation\n");
    return 0;
  }
  
  // Generate testdata
  fprintf(stderr, "C Generate Testdata\n");
  std::vector<PrismCu>  *prisms = generatePrismsFromTestdata(host_mesh_z, host_p_in, host_size_p, host_t_in, host_size_t, host_z_mesh);
  std::vector<PointCu> *samples = generateSamplesFromTestdata(host_mesh_z, host_p_in, host_size_p);
  std::vector<double>    *betas = generateBetasFromTestdata(host_beta_v, host_mesh_z * host_size_t);
  std::vector<double>      *ase = new std::vector<double>(samples->size(), 0);
  rays_total = (unsigned)pow(2,17);

  // Run 
  unsigned i;
 for(i=1; i < argc; ++i){
    if(strncmp(argv[i], "--rays=", 6) == 0){
      const char* pos = strrchr(argv[i],'=');
      rays_total = atoi(pos+1);
    }
  }

  for(i=1; i < argc; ++i){
    if(strncmp(argv[i], "--mode=", 6) == 0){
      if(strstr(argv[i], "bruteforce_gpu") != 0){
  	runtime = runAseBruteforceGpu(samples, prisms, betas, ase, threads, blocks, rays_total);
	strcpy(runmode, "Bruteforce GPU");
	break;
      }
      else if(strstr(argv[i], "naive_ray_propagation") != 0){
	runtime = runNaiveRayPropagation(ase,threads, blocks, rays_total);
	strcpy(runmode, "Naive Ray Propagation GPU");
	break;
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
  fprintf(stderr, "\n");
  fprintf(stderr, "C Statistics\n");
  fprintf(stderr, "C Prism             : %d\n", (int) prisms->size());
  fprintf(stderr, "C Samples           : %d\n", (int) samples->size());
  fprintf(stderr, "C Rays/Sample       : %d\n", rays_total / samples->size());
  fprintf(stderr, "C Rays Total        : %d\n", rays_total);
  fprintf(stderr, "C GPU Blocks        : %d\n", blocks);
  fprintf(stderr, "C GPU Threads/Block : %d\n", threads);
  fprintf(stderr, "C GPU Threads Total : %d\n", threads * blocks);
  fprintf(stderr, "C Runmode           : %s \n", runmode);
  fprintf(stderr, "C Runtime           : %f s\n", runtime / 1000.0);
  fprintf(stderr, "\n");

  return 0;
}


