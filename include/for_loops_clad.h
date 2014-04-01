/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 *
 * This file is part of HASENonGPU
 *
 * HASENonGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HASENonGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HASENonGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */


#include <mesh.h>

#ifndef for_loops_clad_H
#define for_loops_clad_H
float forLoopsClad(
	std::vector<double> *dndtAse,
	unsigned &raysPerSample,
	Mesh *mesh,
	double *betaCells,
	float hostNTot,
	double hostSigmaA,
	double hostSigmaE,
	unsigned hostNumberOfPoints,
	unsigned hostNumberOfTriangles,
	unsigned hostNumberOfLevels,
	float hostThicknessOfPrism,
	float hostCrystalFluorescence	);

#endif /* for_loops_clad_H */
