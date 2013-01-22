__global__ void trace_on_prisms(PrismCu* prisms, const unsigned max_prisms, RayCu* rays, const unsigned max_rays_per_sample, PointCu *samples, const unsigned blocks_per_sample);

__global__ void trace_on_prisms(PrismCu* prisms, const unsigned max_prisms, RayCu* rays, const unsigned max_rays_per_sample, PointCu *samples, const unsigned blocks_per_sample){
  // Cuda ids
  //unsigned tid = threadIdx.x;
  unsigned bid = blockIdx.x + blockIdx.y * gridDim.x;
  unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;

  // Local data
  unsigned prism_i;
  unsigned sample_i = bid / blocks_per_sample;
  RayCu ray = rays[gid];
  unsigned beta_per_ray = 1;
  unsigned importance_per_prism = 1;
  
  // Calculations
  __syncthreads();
  for(prism_i = 0; prism_i < max_prisms; ++prism_i){
    float distance = fabs(collide_prism_gpu(prisms[prism_i], ray));
    ray.P.w += distance * beta_per_ray;
    __syncthreads();	
  }
  __syncthreads();

  // Check Solution
  /* if(fabs(ray.P.w - distance_gpu(ray.P, ray.direction)) > 0.00001){ */
  /*   printf("\033[31;1m[Error]\033[m Sample %d Ray %d with wrong distance real_distance(%f) != sum_distance(%f)\n", sample_i, gid, distance_gpu(ray.P, ray.direction), ray.P.w); */
  /*   return; */
  /* } */

  // Copy data to global
  rays[gid].P.w = ray.P.w;
  //atomicAdd(&(samples[sample_i].w), (ray.P.w * importance_per_prism));

}