#include "datatypes.h"
#ifndef GENERATE_TESTDATA_H
#define GENERATE_TESTDATA_H
std::vector<triangle_cu> generate_triangles(int height, int width, float level);
std::vector<prism_cu> generate_prisms(int height, int width, float level);
std::vector<ray_cu> generate_rays(int height, int width, int level, unsigned max_rays);
std::vector<point_cu> generate_samples(int height, int width, int level);
std::vector<ray_cu> generate_sample_rays(int height, int width, int level, unsigned max_rays, point_cu sample);
ray_cu   generate_ray(int height, int weight, int level);

std::vector<triangle_cu> generate_triangles(int height, int weight, float level){
  int h,w;
  std::vector<triangle_cu> triangles;
  for(h = 0; h < height; ++h){
    for(w = 0; w < weight; ++w){
      triangle_cu t1 = {
	{float(h), float(w), level, 1},
	{float(h), float(w+1), level, 1},
	{float(h+1), float(w), level, 1}};
      triangle_cu t2 = {
	{float(h), float(w+1), level, 1},
	{float(h+1), float(w+1), level, 1},
	{float(h+1), float(w), level, 1}};
      triangles.push_back(t1);
      triangles.push_back(t2);

    }

  }

  return triangles;
}

std::vector<prism_cu> generate_prisms(int height, int width, float level){
  int h,w,l;
  std::vector<prism_cu> prisms;
  for(l = 0; l < level; ++l){
    for(h = 0; h < height; ++h){
      for(w = 0; w < width; ++w){
	triangle_cu a1 = {
	  {float(h), float(w), l, l+1},
	  {float(h), float(w+1), l, l+1},
	  {float(h+1), float(w), l, l+1}};
	triangle_cu b1 = {
	  {float(h), float(w+1), l, 1+1},
	  {float(h+1), float(w+1), l, 1+1},
	  {float(h+1), float(w), l, 1+1}};
      
	prism_cu pr1 = {a1};
	prism_cu pr2 = {b1};

	prisms.push_back(pr1);
	prisms.push_back(pr2);

      }

    }

  }

  return prisms;
}

ray_cu generate_ray(const int heigth, const int width, const int level){
  float rand_heigth = float(rand() % heigth) + (rand() / (float) RAND_MAX);
  float rand_width  = float(rand() % width ) + (rand() / (float) RAND_MAX);
  float rand_level  = float(rand() % level ) + (rand() / (float) RAND_MAX);

  float dir_x = (rand() / (float) RAND_MAX);
  float dir_y = (rand() / (float) RAND_MAX);
  float dir_z = (rand() / (float) RAND_MAX);

  ray_cu r = {
    {rand_heigth, rand_width, rand_level, 0},
    {dir_x, dir_y, dir_z, 0}};
  return r;
}

std::vector<ray_cu> generate_rays(const int height, const int width, const int level, const unsigned max_rays){
  std::vector<ray_cu> rays;
  unsigned ray_i;
  for(ray_i = 0; ray_i < max_rays; ++ray_i){
    ray_cu ray = generate_ray(height, width, level);
    rays.push_back(ray);
  }
  return rays;
}

std::vector<point_cu> generate_samples(int height, int width, int level){
  std::vector<point_cu> sample_points;
  int h,w,l;
  for(l = 0; l <= level; ++l){
    for(h = 0; h <= height; ++h){
      for(w = 0; w <= width; ++w){
	point_cu p = {float(h), float(w), float(l), 0};
	sample_points.push_back(p);
      }
    }
  }
  return sample_points;
}

std::vector<ray_cu> generate_sample_rays(int height, int width, int level, unsigned max_rays, point_cu sample){
  std::vector<ray_cu> rays;
  unsigned ray_i;
  for(ray_i = 0; ray_i < max_rays; ++ray_i){
    ray_cu ray = generate_ray(height, width, level);
    ray.direction.x  =  ray.P.x;
    ray.direction.y  =  ray.P.y;
    ray.direction.z  =  ray.P.z;
    ray.P = sample;
    rays.push_back(ray);
  }
  return rays;
}

#endif /* GENERATE_TESTDATA_H */
