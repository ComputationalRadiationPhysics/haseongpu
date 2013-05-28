#include "datatypes.h"

#ifndef GENERATE_TESTDATA_H
#define GENERATE_TESTDATA_H
std::vector<PointCu> *generateSamplesFromTestdata(unsigned levels, std::vector<double>* points, unsigned numPoints);
std::vector<PrismCu> *generatePrismsFromTestdata(unsigned levels, std::vector<double>* points, unsigned numPoints, std::vector<unsigned> *triangles, unsigned numTriangles, double thickness);
std::vector<double> *generateBetasFromTestdata(std::vector<double>* betaValues, unsigned numBetaValue);

#endif /* GENERATE_TESTDATA_H */
