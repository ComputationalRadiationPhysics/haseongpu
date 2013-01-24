#include "datatypes.h"
#include <vector>

#ifndef GENERATE_TESTDATA_H
#define GENERATE_TESTDATA_H
std::vector<TriangleCu> generate_triangles(int length, int width, float level); // not used
std::vector<PrismCu> *generate_prisms(int length, int width, float level);
std::vector<RayCu> generate_rays(int length, int width, int level, unsigned max_rays); // not used
std::vector<PointCu> *generate_samples(int length, int width, int level);
std::vector<RayCu> *generate_sample_rays(int length, int width, int level, unsigned max_rays, std::vector<PointCu> samples);
RayCu generate_ray(int length, int weight, int level);

std::vector<TriangleCu> generate_triangles(int length, int weight, float level){
  int h,w;
  std::vector<TriangleCu> triangles;
  for(h = 0; h < length; ++h){
    for(w = 0; w < weight; ++w){
      TriangleCu t1 = {
	{float(h), float(w), level, 1},
	{float(h), float(w+1), level, 1},
	{float(h+1), float(w), level, 1}};
      TriangleCu t2 = {
	{float(h), float(w+1), level, 1},
	{float(h+1), float(w+1), level, 1},
	{float(h+1), float(w), level, 1}};
      triangles.push_back(t1);
      triangles.push_back(t2);

    }

  }

  return triangles;
}

std::vector<PrismCu> *generate_prisms(int length, int width, float level){
  int h,w,l;
  std::vector<PrismCu> *prisms = new std::vector<PrismCu>;
  for(l = 0; l < level; ++l){
    for(h = 0; h < length; ++h){
      for(w = 0; w < width; ++w){
	TriangleCu a1 = {
	  {float(h), float(w), l, 1},
	  {float(h), float(w+1), l, 1},
	  {float(h+1), float(w), l, 1}};
	TriangleCu b1 = {
	  {float(h), float(w+1), l, 1},
	  {float(h+1), float(w+1), l, 1},
	  {float(h+1), float(w), l, 1}};
      
	PrismCu pr1 = {a1};
	PrismCu pr2 = {b1};

	prisms->push_back(pr1);
	prisms->push_back(pr2);

      }

    }

  }

  return prisms;
}

std::vector<PrismCu> *generatePrismsFromTestdata(unsigned levels, double *points, unsigned *triangles, numTriangles, unsigned thickness){
  std::vector<PrismCu> *prisms = new std::vector<PrismCu>;
  unsigned triangle_i;
  unsigned level_i;

  for(level_i = 0; level_i < levels; ++level_i){
  for(triangle_i = 0; triangle_i < numTriangles; ++triangle_i){
    unsigned a = triangles[triangle_i * 3];
    unsigned b = triangles[triangle_i * 3 + 1];
    unsigned c = triangles[triangle_i * 3 + 2];

    PointCu A = {
      points[a * 2],
      points[a * 2 + 1],
      level_i,
      thickness};

    PointCu B = {
      points[b * 2],
      points[b * 2 + 1],
      level_i,
      thickness};

    PointCu C = {
      points[c * 2],
      points[c * 2 + 1],
      level_i,
      thickness};

    TriangleCu triangle = {A, B, C};
    PrismCu prism = {triangle};
    prisms->push_back(prism);
  }
  }

  return prisms;
}

RayCu generate_ray(const int heigth, const int width, const int level){
  float rand_heigth = float(rand() % heigth) + (rand() / (float) RAND_MAX);
  float rand_width  = float(rand() % width ) + (rand() / (float) RAND_MAX);
  float rand_level  = float(rand() % level ) + (rand() / (float) RAND_MAX);

  float dir_x = (rand() / (float) RAND_MAX);
  float dir_y = (rand() / (float) RAND_MAX);
  float dir_z = (rand() / (float) RAND_MAX);

  RayCu r = {
    {rand_heigth, rand_width, rand_level, 0},
    {dir_x, dir_y, dir_z, 0}};
  return r;
}

std::vector<RayCu> generate_rays(const int length, const int width, const int level, const unsigned max_rays){
  std::vector<RayCu> rays;
  unsigned ray_i;
  for(ray_i = 0; ray_i < max_rays; ++ray_i){
    RayCu ray = generate_ray(length, width, level);
    rays.push_back(ray);
  }
  return rays;
}

std::vector<PointCu> *generate_samples(int length, int width, int level){
  std::vector<PointCu> *sample_points = new std::vector<PointCu>;
  int h,w,l;
  for(l = 0; l <= level; ++l){
    for(h = 0; h <= length; ++h){
      for(w = 0; w <= width; ++w){
	PointCu p = {float(h), float(w), float(l), 0};
	sample_points->push_back(p);
      }
    }
  }
  return sample_points;
}

std::vector<RayCu> *generate_sample_rays(int length, int width, int level, unsigned max_rays, std::vector<PointCu> *samples){
  std::vector<RayCu> *rays = new std::vector<RayCu>;
  unsigned ray_i, sample_i;
  for(sample_i = 0; sample_i < samples->size(); ++sample_i){
    for(ray_i = 0; ray_i < max_rays; ++ray_i){
      RayCu ray = generate_ray(length, width, level);
      ray.direction.x  =  ray.P.x;
      ray.direction.y  =  ray.P.y;
      ray.direction.z  =  ray.P.z;
      ray.P = samples->at(sample_i);
      ray.P.w = 0.f;
      rays->push_back(ray);
    }

  }
  return rays;
}



#endif /* GENERATE_TESTDATA_H */
