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
  const unsigned max_rays = 1024;
  const unsigned max_triangles = 8;
  const unsigned depth  = 2;
  const unsigned length = ceil(sqrt(max_triangles / 2));
  const int threads = 256;
  unsigned ray_i, prism_i, sample_i, i;
  float runtime_gpu = 0.0;
  float runtime_cpu = 0.0;
  float runtime = 0.0;
  bool use_cpu = true;
  bool use_gpu = true;
  
  // Parse Commandline
  if(argc <= 1){
    fprintf(stderr, "C No commandline arguments found\n");
    fprintf(stderr, "C Usage    : ./octrace --mode=[runmode]\n");
    fprintf(stderr, "C Runmodes : bruteforce_cpu\n");
    fprintf(stderr, "             bruteforce_gpu\n");
    return 0;
  }

  /* for(i=1; i < argc; ++i){ */
  /*   if(strncmp(argv[i], "--mode=", 6) == 0){ */
  /*     if(strstr(argv[i], "bruteforce_cpu") != 0){ */
  /* 	runtime_cpu = run_ase_bruteforce_cpu(samples, prisms, rays, ase_cpu); */
  /*     }else */
  /*     if(strstr(argv[i], "bruteforce_gpu") != 0){ */
  /* 	runtime_gpu = run_ase_bruteforce_gpu(samples, prisms, rays, ase_gpu, threads); */
  /*     } */


  /*   } */
       
  /* } */
  
  // Generate testdata
  fprintf(stderr, "C Generate Testdata\n");
  std::vector<PrismCu> *prisms  = generate_prisms(length, length, depth);
  std::vector<PointCu> *samples = generate_samples(length, length, depth);
  std::vector<RayCu> *rays      = generate_sample_rays(length, length, depth, max_rays, samples);
  std::vector<float> *ase_cpu   = new std::vector<float>(samples->size(), 0);
  std::vector<float> *ase_gpu   = new std::vector<float>(samples->size(), 0);
  
  // CPU Raytracing
  if(use_cpu){
    runtime_cpu = run_ase_bruteforce_cpu(samples, prisms, rays, ase_cpu);
  }
   
  // GPU Raytracing
  if(use_gpu){
    runtime_gpu = run_ase_bruteforce_gpu(samples, prisms, rays, ase_gpu, threads);
  }

  // Evaluate device data
  if(use_gpu && use_cpu){
    for(sample_i = 0; sample_i < samples->size(); ++sample_i){
      if(fabs(ase_cpu->at(sample_i) - ase_gpu->at(sample_i)) < 0.1){
  	fprintf(stderr, "CPU == GPU: Sample %d with value %f \n", sample_i, ase_cpu->at(sample_i));

      }
      else{
  	fprintf(stderr, "\033[31;1m[Error]\033[m CPU(%f) != GPU(%f) on sample %d\n", ase_cpu->at(sample_i), ase_gpu->at(sample_i), sample_i);

      }

    }

    /* for(ray_i = 0; ray_i < rays.size(); ++ray_i){ */
    /*   if(fabs(ray_data[ray_i] - h_rays[ray_i].P.w) < 0.00001){ */
    /* 	//fprintf(stderr, "CPU == GPU: Ray %d with value %f \n", ray_i, ray_data[ray_i]); */

    /*   } */
    /*   else{ */
    /* 	fprintf(stderr, "\033[31;1m[Error]\033[m CPU(%f) != GPU(%f) on ray %d\n", ray_data[ray_i], h_rays[ray_i].P.w, ray_i); */

    /*   } */
    /* } */
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
  fprintf(stderr, "C Runtime_GPU       : %f s\n", runtime_gpu / 1000.0);
  fprintf(stderr, "C Runtime_CPU       : %f s\n", runtime_cpu / 1000.0);
  fprintf(stderr, "C Speedup CPU/GPU   : %.1f\n", runtime_cpu / runtime_gpu);
  fprintf(stderr, "\n");





  // Cleanup
  /* cudaFreeHost(h_rays); */
  /* cudaFreeHost(h_prisms); */
  /* cudaFreeHost(h_samples); */
 
  return 0;
}

