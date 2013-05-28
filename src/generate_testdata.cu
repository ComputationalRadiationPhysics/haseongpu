#include "datatypes.h"
#include <vector> /* vector */
#include "assert.h" /* assert */
#include "stdio.h" /* printf */

std::vector<PointCu> *generateSamplesFromTestdata(unsigned levels, std::vector<double>* points, unsigned numPoints){
  std::vector<PointCu> *p = new std::vector<PointCu>;
  unsigned level_i;
  unsigned point_i; 

  for(level_i = 0; level_i < levels; ++level_i){
    for(point_i = 0; point_i < numPoints; ++point_i){
      PointCu point = {
  	points->at(point_i),
  	points->at(numPoints + point_i),
  	level_i,
  	0};
      p->push_back(point);
    }

  }
  return p;
}

std::vector<PrismCu> *generatePrismsFromTestdata(unsigned levels, std::vector<double>* points, unsigned numPoints, std::vector<unsigned> *triangles, unsigned numTriangles, double thickness){
  std::vector<PrismCu> *prisms = new std::vector<PrismCu>;
  unsigned triangle_i;
  unsigned level_i;
  for(level_i = 0; level_i < levels-1; ++level_i){
    for(triangle_i = 0; triangle_i < numTriangles; ++triangle_i){
      unsigned a = triangles->at(triangle_i);
      unsigned b = triangles->at(numTriangles + triangle_i);
      unsigned c = triangles->at(2 * numTriangles + triangle_i);

      PointCu A = {
	points->at(a),
	points->at(numPoints + a),
	level_i,
	thickness};

      PointCu B = {
	points->at(b),
	points->at(numPoints + b),
	level_i,
	thickness};

      PointCu C = {
	points->at(c),
	points->at(numPoints + c),
	level_i,
	thickness};

      TriangleCu triangle = {A, B, C};
      PrismCu prism = {triangle};
      prisms->push_back(prism);
    }
  }
  
  return prisms;
}

std::vector<double> *generateBetasFromTestdata(std::vector<double>* betaValues, unsigned numBetaValue){
  std::vector<double> *betas = new std::vector<double>;
  
  unsigned beta_i;
  for(beta_i = 0; beta_i < numBetaValue; ++beta_i){
    betas->push_back(betaValues->at(beta_i));
  }

  return betas;

}
