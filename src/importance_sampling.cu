#include <mesh.h>
#include <stdio.h>
#include <propagate_ray.h>
#include <geometry.h>

// ##############################################################
// # Reconstruction                                             #
// ##############################################################
void importanceSamplingNew(Point samplePoint, Mesh mesh, unsigned raysPerSample, double sigmaA, double sigmaE, double nTot,  double *importance, unsigned *raysPerPrism){
  Triangle *triangles = mesh.triangles;
  unsigned numberOfLevels = mesh.numberOfLevels;
  unsigned numberOfTriangles = mesh.numberOfTriangles;
  float thickness = mesh.thickness;

  int raysLeft = 0;
  int raysDump = 0;
  double sumPhi = 0;
  double surfaceTotal = 0;
  double gain = 0;
  Ray ray;
  Point startPoint;
  Triangle startTriangle;

  // Calculate importance by propagation from trianglecenter to every other center
  for(unsigned triangle_i = 0; triangle_i < numberOfTriangles; ++triangle_i){
    for(unsigned level_i = 0; level_i < numberOfLevels - 1; ++level_i){
      startTriangle = triangles[triangle_i];
      startPoint.x = startTriangle.center.x;
      startPoint.y = startTriangle.center.y;
      startPoint.z = (level_i + 0.5) * thickness;
      // DEBUG
      // printf("\nstartpoint x %f\n", startPoint.x );
      // printf("startpoint y %f\n",  startPoint.y);
      // printf("startpoint z %f\n",  startPoint.z);
      // printf("endpoint x %f\n", samplePoint.x);
      // printf("endpoint y %f\n",  samplePoint.y);
      // printf("endpoint z %f\n", samplePoint.z);
      ray = generateRay(startPoint, samplePoint);

      gain = propagateRay(ray, level_i, startTriangle, triangles, sigmaA, sigmaE, nTot, thickness);

      importance[triangle_i + level_i * numberOfTriangles] = startTriangle.betaValues[level_i] * gain;
      sumPhi += importance[triangle_i + level_i * numberOfTriangles];

    }
    surfaceTotal += triangles[triangle_i].surface;
  }

  // Calculate number of rays/prism
  for(unsigned triangle_i = 0; triangle_i < numberOfTriangles; ++triangle_i){
    for(unsigned level_i = 0; level_i < numberOfLevels - 1; ++level_i){
      raysPerPrism[triangle_i + level_i * numberOfTriangles] =  (unsigned)(floor(importance[triangle_i + level_i * numberOfTriangles] / sumPhi * raysPerSample));
      raysDump +=  raysPerPrism[triangle_i + level_i * numberOfTriangles];
      raysDump +=  raysPerPrism[0];
    }

  }

  raysLeft = raysPerSample - raysDump;

  // TODO What happens with random failure ?
  // TODO Distribute the remaining rays randomly
  for (int i_r=0; i_r < raysLeft; i_r++){
    int rand_t = (int )(rand() % numberOfTriangles);
    int rand_z = (int )(rand() % (numberOfLevels-1));
    raysPerPrism[rand_t + rand_z * numberOfTriangles]++;

  }


  //  Now think about the mount of rays which would come out of this volume(surface)
  //  dividing this number with the new amount of rays gives the final importance weight for this area!
  // for (int triangle_i=0; triangle_i < numberOfTriangles; ++triangle_i){
  //   for (int level_i=0; level_i < numberOfLevels; ++level_i){
  //     if (raysPerPrism[triangle_i + (level_i * numberOfTriangles)] > 0){
  // 	importance[triangle_i + (level_i * numberOfTriangles)] = raysPerSample * triangles[triangle_i].surface / surfaceTotal / raysPerPrism[triangle_i + (level_i * numberOfTriangles)];

  //     }
  //     else{
  // 	importance[triangle_i + (level_i * numberOfTriangles)] = 0; 

  //     }

  //   }

  // }

}


// #################################################
// # Old Code                                      #
// #################################################

double propagateRayHost (double x_pos, 
			 double y_pos, 
			 double z_pos, 
			 double x_dest, 
			 double y_dest, 
			 double z_dest, 
			 int t_start, 
			 int mesh_start, 
			 double *p_in,
			 double *n_x,
			 double *n_y,
			 int *n_p,
			 int *neighbors,
			 int *forbidden,
			 double *beta_v,
			 unsigned numberOfPoints,
			 unsigned numberOfTriangles,
			 float thicknessOfPrism,
			 float sigmaA,
			 float sigmaE,
			 float nTot
			 ){
  //    in first try no reflections
  //    calculate the vector and make the the calculation, which surface would be the shortest to reach
  //    then get the length, make the integration, get the information about the next cell out of the array
  //    set the point to the surface (this surface is "forbidden" in the calculations)
  //    proceed until you hit a the point or the surface
  //    if you are closer then "small" stop and return the value
  double vec_x, vec_y,vec_z, norm;
  double distance, length, length_help, distance_total;
  double gain=1;
  double nominator, denominator;
  int tri, cell_z; // the current triangle number and position concerning the z's
  int decider; // which one is the shortest - info
  int tri_next, cell_z_next, forb, forb_dump;
  unsigned size_p = numberOfPoints;
  unsigned N_cells = numberOfTriangles;
  float z_mesh = thicknessOfPrism;
  float sigma_a = sigmaA;
  float sigma_e = sigmaE;
  float N_tot = nTot;

	
  //    initial positions
  tri = t_start;
  cell_z = mesh_start;
    
    
  //    definition of the vectors without reflections
  vec_x = (x_dest - x_pos);
  vec_y = (y_dest - y_pos);
  vec_z = (z_dest - z_pos);
    
  norm = sqrt(vec_x*vec_x+vec_y*vec_y+vec_z*vec_z);
    
  vec_x = vec_x/norm;
  vec_y = vec_y/norm;
  vec_z = vec_z/norm;
    
  //    now calculate the length to travel
  distance = sqrt((x_dest - x_pos)*(x_dest - x_pos)+(y_dest - y_pos)*(y_dest - y_pos)+(z_dest - z_pos)*(z_dest - z_pos));
  distance_total = distance;
  // does this make sense?
  length = distance;
    
  forb = -1;
	
  //	mexPrintf("Propagation called");
  //    mexEvalString("drawnow;");
    
  //    the ray has to be set to be ALIVE before!
  //    now do the unlimited for loop - break!!!
  for(;;)
    {
	
      //	  mexPrintf("Propagation for part called\n\n");
      //    mexEvalString("drawnow;");
      //        definition for decider
      //        0,1,2: int for the neighbors
      //        3: hor plane up
      //        4: hor plane down
        
      //        at first set the decider = -1;
      decider = -1;
      length = distance;
		
		
      //        mexPrintf("forb: %i\n",forb);
      //        mexEvalString("drawnow;");
        
      //        try the triangle faces
      //        remember the correlation between the normals and the points
      //        n1: p1-2, n2: p1-3, n3:p2-3
      //        the third coordinate (z) of the particpating points for the surfaces can be set to be z=0, 
      //        as everything uses triangular "tubes/prisms", as well as n_z=0 in this case!
      if (forb != 0){
	nominator = (n_x[tri]*p_in[n_p[tri]] + n_y[tri]*p_in[n_p[tri]+size_p]) - (n_x[tri]*x_pos + n_y[tri]*y_pos);
	denominator = n_x[tri]*vec_x + n_y[tri]*vec_y;
	if (denominator != 0.0)
	  {
	    length_help = nominator/denominator;
	    if (length_help < length && length_help > VERY_SMALL)
	      {
		length = length_help;
		decider = 0;
		forb_dump = (forbidden[tri]);
	      }
	  }
      }
        
      if (forb != 1){
	nominator = (n_x[tri+N_cells]*p_in[n_p[tri+N_cells]] + n_y[tri+N_cells]*p_in[n_p[tri+N_cells]+size_p]) - (n_x[tri+N_cells]*x_pos + n_y[tri+N_cells]*y_pos);
	denominator = n_x[tri+N_cells]*vec_x + n_y[tri+N_cells]*vec_y;
	if (denominator != 0.0)
	  {
	    length_help = nominator/denominator;
	    if (length_help < length && length_help > VERY_SMALL)
	      {
		length = length_help;
		decider = 1;
		forb_dump = (forbidden[tri+N_cells]);
	      }
	  }
      }
        
      if (forb !=2){
	nominator = (n_x[tri+2*N_cells]*p_in[n_p[tri+2*N_cells]] + n_y[tri+2*N_cells]*p_in[n_p[tri+2*N_cells]+size_p]) - (n_x[tri+2*N_cells]*x_pos + n_y[tri+2*N_cells]*y_pos);
	denominator = n_x[tri+2*N_cells]*vec_x + n_y[tri+2*N_cells]*vec_y;
	if (denominator != 0.0)
	  {
	    length_help = nominator/denominator;
	    if (length_help < length && length_help > VERY_SMALL)
	      {
		length = length_help;
		decider = 2;
		forb_dump = (forbidden[tri+2*N_cells]);
	      }
	  }
      }
        
      //        try the horizontal planes, which one is the shortest, n_x and n_y are zero!, n_z =1!
      //        at first the upper plane
      if (forb != 3){
	nominator = (cell_z+1)*z_mesh - z_pos;
	denominator = z_pos*vec_z;
	if (denominator != 0.0)
	  {
	    length_help = nominator/denominator;
	    if (length_help < length && length_help > VERY_SMALL)
	      {
		length = length_help;
		decider = 3;
		forb_dump = 4; // you are not allowed to go down in the next step
	      }
	  }
      }
        
      //        next is the lower plane
      if (forb != 4){
	nominator = (cell_z)*z_mesh - z_pos;
	denominator = z_pos*vec_z;
            
	if (denominator != 0.0)
	  {
	    length_help = nominator/denominator;
	    if (length_help < length && length_help > VERY_SMALL)
	      {
		length = length_help;
		decider = 4;
		forb_dump = 3; // you are not allowed to go up in the next step
	      }
	  }
      }
        
      forb = forb_dump;
		
        
      //        now make a switch to differ the different cases
      switch(decider){
                
      case 0:
	//                this is the case for the intersection with the first choice triangle-surface
	tri_next = neighbors[tri];
	cell_z_next = cell_z;
	break;
                
      case 1:
	//                second triangle surface
	tri_next = neighbors[tri+N_cells];
	cell_z_next = cell_z;
	break;
                
      case 2:
	//                third triangle surface
	tri_next = neighbors[tri+2*N_cells];
	cell_z_next = cell_z;
	break;
                
      case 3:
	//                go one plane up
	tri_next = tri;
	cell_z_next = cell_z + 1;
	break;
                
      case 4:
	//                go one plane down
	tri_next = tri;
	cell_z_next = cell_z - 1;
	break;
                
      default:
	//                make an error statement
	break;
      }
        
      //        now we know where to go, let's make the integration
      //        take the beta_v[tri+cell_z*N_cells] 

      //		  at this position do the decision whether it is a gain part or cladding
      //		  it might be absorbing or amplifying, for the cladding only absorbing
      //		  a simple "if then"

      gain = gain * exp(N_tot*(beta_v[tri+cell_z*N_cells]*(sigma_e + sigma_a)-sigma_a)*length);

      //        gain = LineIntegralMCRK4_S(3, tri, cell_z, gain, length);
        
      //        after integration make the propagation
        
      //        mexPrintf("Distance: %f, Length: %f\n",distance, length);
      //        mexPrintf("decider: %i, forbidden: %i\n",decider, forb);
      //        mexPrintf("vec_x: %f, vec_y: %f, vec_z: %f\n", vec_x, vec_y, vec_z);
      //        mexPrintf("current_x: %f current_y: %f current_z: %f\n", x_pos, y_pos, z_pos);
      //        mexPrintf("tri: %i, tri_next: %i, cell_z: %i, cell_next: %i\n", tri, tri_next, cell_z, cell_z_next);
      //        mexEvalString("drawnow;");
      //        str=mxCreateString("Press a key");
      //        mexCallMATLAB(1,&dump,1,&str,"input"); 
      //        str and dump should be defined to be a *mxArray and don't forget to kill them at the end
        
      distance -= length;
        
      //        return 1;
      //        
        
      x_pos = x_pos + length*vec_x;
      y_pos = y_pos + length*vec_y;
      z_pos = z_pos + length*vec_z;
        
      if (abs(distance)< SMALL)
        {
	  break;
        }
        
        
      //        now set the next cell
      tri = tri_next;
      cell_z = cell_z_next;      
        
      //        break;
      //        now we should make the integration routine
    }
    
  gain /= (distance_total*distance_total);

  return gain;
}

/**
 * calculate the gain from the centers of each of the boxes to the observed point
 * calculate the gain and make a "mapping"
 * receipt: pick the point in the center of one cell, 
 * calculate the gain from this point to the observed point,
 * estimate the inner part of the Phi_ASE - Integral,
 * scale the amount of rays proportionally with it
 * sum the amount of rays and scale it to Int=1, which gives the inverse weights
 * the number of rays is determined via floor(), with ceil(), zero-redions could be added
 * use the routine "propagation"!, test: no reflections, just exponential
 *
 **/
void importanceSampling(int point,
	     int startLevel,
	     double *importance,
	     unsigned *numberOfImportantRays,
	     double *points,
	     double *xOfNormals,
	     double *yOfNormals,
	     int *positionsOfNormalVectors,
	     int *neighbors,
	     int *forbidden,
	     double *betaValues,
	     double *xOfTriangleCenter,
	     double *yOfTriangleCenter,
	     float *surface,
	     unsigned raysPerSample,
	     unsigned numberOfPoints,
	     unsigned numberOfLevels,
	     unsigned numberOfTriangles,
	     float thicknessOfPrism,
	     float sigmaA,
	     float sigmaE,
	     float nTot
	     )
{
  int raysLeft;
  int raysDump;
  double sumPhi;
  double surfaceTotal;
  double xPos, yPos, zPos;
  double prop;

  raysDump = 0;
  sumPhi = 0.0;
  surfaceTotal = 0.0;
  xPos = points[point];
  yPos = points[point + numberOfPoints];
  zPos = startLevel * thicknessOfPrism;

  // Calculate importance by propagation from trianglecenter to every other center
  for (int i_t=0; i_t < numberOfTriangles; ++i_t){
    for (int i_z=0; i_z < (numberOfLevels-1); ++i_z){
      prop = propagateRayHost(xOfTriangleCenter[i_t], yOfTriangleCenter[i_t], 
      			       thicknessOfPrism * (i_z+0.5),  xPos, yPos, zPos, i_t, i_z, 
      			       points, xOfNormals, yOfNormals, positionsOfNormalVectors, 
      			       neighbors, forbidden, betaValues,
      			       numberOfPoints, numberOfTriangles, thicknessOfPrism,
      			       sigmaA, sigmaE, nTot
      			       );

      importance[i_t + i_z * numberOfTriangles] = betaValues[i_t + i_z * numberOfTriangles]*(prop);
      sumPhi += importance[i_t + i_z * numberOfTriangles];

    }
    surfaceTotal += surface[i_t];

  }

  // Calculate number of rays/prism
  for (int i_t=0; i_t < numberOfTriangles; ++i_t){
    for (int i_z=0; i_z < (numberOfLevels-1); ++i_z){
      numberOfImportantRays[i_t + i_z * numberOfTriangles] = (unsigned)(floor(importance[i_t + i_z * numberOfTriangles] / sumPhi * raysPerSample));
      raysDump +=  numberOfImportantRays[i_t + i_z * numberOfTriangles];
      //fprintf(stderr, "[%d][%d] i: %.20f n: %d\n", i_z, i_t, importance[i_t + i_z*numberOfTriangles], numberOfImportantRays[i_t + i_z*numberOfTriangles]);
    }

  }
  raysLeft = raysPerSample - raysDump;

  // TODO What happens with random failure ?
  // Distribute the remaining rays randomly
  for (int i_r=0; i_r < raysLeft; i_r++){
    int rand_t = (int )(rand() % numberOfTriangles);
    int rand_z = (int )(rand() % (numberOfLevels-1));
    numberOfImportantRays[rand_t + rand_z * numberOfTriangles]++;

  }

  //  Now think about the mount of rays which would come out of this volume(surface)
  //  dividing this number with the new amount of rays gives the final importance weight for this area!
  for (int i_t=0; i_t<numberOfTriangles; ++i_t){
    for (int i_z=0; i_z<(numberOfLevels-1); ++i_z){
      if (numberOfImportantRays[i_t + i_z*numberOfTriangles] > 0){
  	importance[i_t + i_z*numberOfTriangles] = raysPerSample * surface[i_t] / surfaceTotal / numberOfImportantRays[i_t + i_z*numberOfTriangles];

      }
      else{
  	importance[i_t + i_z*numberOfTriangles] = 0; 

      }


    }

  }


}
