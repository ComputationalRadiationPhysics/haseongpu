// Libraries
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "vector_types.h"
#include "assert.h"
#include "string.h"
#include <vector>
#include "curand_kernel.h"
#include "datatypes.h"

// User header files
#include "geometry.h"
#include "datatypes.h"
#include "generate_testdata.h"
#include "print.h"
#include "geometry_gpu.h"
#include "ase_bruteforce_kernel.h"
#include "ase_bruteforce_cpu.h"
#include "testdata.h"

int main(int argc, char **argv){
  const unsigned rays_per_sample = 100000;
  const unsigned max_triangles = 2;
  const unsigned depth  = 2;
  const unsigned length = ceil(sqrt(max_triangles / 2));
  const int threads = 256;
  char runmode[20];
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
  //std::vector<PrismCu> *prisms  = generate_prisms(length, length, depth);
  std::vector<PrismCu> *prisms  = generatePrismsFromTestdata(host_z_mesh, host_p_in, host_t_in, host_number_of_triangles, host_mesh_z);
  //std::vector<PointCu> *samples = generate_samples(length, length, depth);
  std::vector<PointCu> *samples = generateSampesFromTestdata(host_z_mesh, host_p_in, host_number_of_points);
  std::vector<RayCu> *rays      = generate_sample_rays(length, length, depth, rays_per_sample, samples);
  std::vector<float> *ase       = new std::vector<float>(samples->size(), 0);

  // Run 
  unsigned i;
  for(i=1; i < argc; ++i){
    if(strncmp(argv[i], "--mode=", 6) == 0){
      if(strstr(argv[i], "bruteforce_cpu") != 0){
  	runtime = runAseBruteforceCpu(samples, prisms, rays, ase);
	strcpy(runmode, "Bruteforce CPU");

      }
      else if(strstr(argv[i], "bruteforce_gpu") != 0){
  	runtime = runAseBruteforceGpu(samples, prisms, rays, ase, threads);
	strcpy(runmode, "Bruteforce GPU");

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
  for(sample_i = 0; sample_i < samples->size(); ++sample_i){
    fprintf(stderr, "C ASE PHI of sample %d: %f\n", sample_i, ase->at(sample_i));

  }

  // Print statistics
  unsigned blocksPerSample = ceil(rays->size() / (threads * samples->size()));
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


