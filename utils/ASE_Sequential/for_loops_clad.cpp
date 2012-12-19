#include <math.h>
#include "mt19937ar.h"
#include <matrix.h>
#include <mex.h>

#define ALIVE 1
#define DEAD 0
#define SMALL 1E-06 // 1µm is considered to be small
#define MAXREFL 5 // number of max amount of reflections <= this has to be done in a better way later (distance defines the max number)

// change 2011/02/23
// added external input of emission and absorption cross section

// change 2011/04/18
// cladding information and calculation added

// change 2011/04/21
// compiler note: mex for_loops_clad.cpp mt19937ar.cpp

//global variables
double *p_in, *beta_p, *beta_v, *n_x, *n_y, *surface, *center_x, *center_y;
int *t_in, *n_p, *neighbors, *forbidden;
int size_p, N_cells, mesh_z, rays;
int NumRays;
double z_mesh=0.0356;
double clad_abs;
int clad_num, *cell_type;


//ray_life determines, if the ray is DEAD or ALIVE
int ray_life;

//  spectral properties
double beta = 0.12;
double sigma_a = 1.16e-21;
double sigma_e = 2.48e-20;
double doping = 5;
double N_tot = 1.38e20;


// function prototypes
double LineIntegralMCRK4_S(int nSteps, int i_t, int i_z, double PhotonsIn, double length);
double propagation(double x_pos, double y_pos, double z_pos, double x_dest, double y_dest, double z_dest, int t_start, int mesh_start, int N_refl);
void importf(int point, int mesh_start, double *importance, int *N_rays, int N_reflections);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double u=0.0, v=0.0, w=0.0;  
  int i=0, j=0, k=0, l=0, iz=0;  
  int t_1,t_2,t_3;
  double p_cx,p_cy, p_cz;
  double x_rand, y_rand, z_rand;
//  double distance;
  double gain;
  
//  int N_rays_dump;
  
//  mex inputs and outputs
//  double *p_in, *beta_p, *beta_v, *n_x, *n_y, *surface, *center_x, *center_y;
//  int *t_in, *neighbors, *n_p;
//  defined global ^^
  mxArray *rand_out, *phi_ASE, *import, *N_r;
  double *rand, *phi, *importance;
  int size_t;
  int *N_rays, *Rays;
  double *dop, *zm;
  double *sa, *se; // doping information input
  double *cladabs; // cladding information
  int *cladnum; //cladding information
  int N_refl; // number of reflections that the beam will do
  

//  be very careful how the data is sent to this file and how it has to be cast!
//  remember, some are int32 and some are double32, some are arrays, some not!
  
// ************* INPUT ORDER BEGIN *************
//  for further informations take a look on "mesh_cyl_rect.m"
  
  p_in = (double *)mxGetData(prhs[0]); //point coordinates in 2D
  t_in = (int *)mxGetData(prhs[1]);  //association triangle-points - c-indexing-sytle!
  beta_p = (double *)mxGetData(prhs[2]); //beta in the points
  beta_v = (double *)mxGetData(prhs[3]); //interpolated volume beta
  n_x = (double *)mxGetData(prhs[4]); //normals of the facets, x-components
  n_y = (double *)mxGetData(prhs[5]); //normals of the facets, y-components
//  n_z = (double *)mxGetData(prhs[X]); //would be facets z-normals, but it exists just in the tetrahedron triangulation
  neighbors = (int *)mxGetData(prhs[6]); //for each cell in t_in, its neighboring cell in plane geometry  - c-indexing-sytle!
  surface = (double *)mxGetData(prhs[7]);  //information about the surfaces of the cells
  center_x = (double *)mxGetData(prhs[8]); //center positions of the cells, x-coordinates
  center_y = (double *)mxGetData(prhs[9]); //center positions of the cells, y_coordinates
  n_p = (int *)mxGetData(prhs[10]); // gives the Index to one of the points in p_in which is in the plane of the normals - c-indexing-sytle!
  forbidden = (int *)mxGetData(prhs[11]);//are the correspondance of the face used in the previous triangle, which must not be tested (rounding errors)
  Rays = (int *)mxGetData(prhs[12]);// amount of test rays per cell
  dop = (double *)mxGetData(prhs[13]);//doping informations
  zm = (double *)mxGetData(prhs[14]);//z_mesh info
  se = (double *)mxGetData(prhs[15]);//emission cross section
  sa = (double *)mxGetData(prhs[16]);//absorption cross section
  cell_type = (int *)mxGetData(prhs[17]);//cell type array with one of them as cladding decided by clad_num
  cladnum = (int *)mxGetData(prhs[18]);//number, which one is the cladding
  cladabs = (double *)mxGetData(prhs[19]);//amount of absorption in the cladding 
  
// *************  INPUT ORDER END  *************
  
//  associate inputs
  size_p = (int )mxGetM(prhs[0]); //number of points
  size_t = (int )mxGetM(prhs[1]); //number of triangles per sheet
  mesh_z = (int )mxGetN(prhs[2]); //number of meshing in z-direction
  
//  mexPrintf("mesh_z: %i\n", mesh_z);
//  mexEvalString("drawnow;");
  
  
  int dim[]={size_t,mesh_z-1};
  
// ************* OUTPUT ORDER BEGIN *************
  rand_out = plhs[0] = mxCreateDoubleMatrix(rays,2,mxREAL);
  phi_ASE = plhs[1] = mxCreateDoubleMatrix(size_p,mesh_z,mxREAL);
  import = plhs[2] = mxCreateDoubleMatrix(size_t,mesh_z-1,mxREAL);
  N_r = plhs[3] = mxCreateNumericArray(2,dim,mxINT32_CLASS,mxREAL);
// *************  OUTPUT ORDER END  *************
  
  
//  associate outputs
  rand = (double *)mxGetData(plhs[0]);
  phi = (double *)mxGetData(plhs[1]);
  importance = (double *)mxGetData(plhs[2]);
  N_rays = (int *)mxGetData(plhs[3]);
  
  NumRays = (int)(Rays[0]);
  N_tot = (double)(dop[0]);
  z_mesh = (double)(zm[0]);
  sigma_a = (double)(sa[0]);
  sigma_e = (double)(se[0]);
  clad_num = (int)(cladnum[0]);
  clad_abs = (double)(cladabs[0]);
  
  N_cells = size_t;
  mexPrintf("NumRays: %i\n", NumRays);
  mexEvalString("drawnow;");
   
//  print the first 3 values
//  mexPrintf("the first point is: (%1.5f, %1.5f)\n",p_in[0],p_in[size_p]);
//  mexPrintf("the second point is: (%1.5f, %1.5f)\n",p_in[1],p_in[size_p+1]);
//  mexPrintf("the third point is: (%1.5f, %1.5f)\n\n",p_in[2],p_in[size_p+2]);
//  
  
//  make drawnow to force MatLab to show the bash output during run
//  mexEvalString("drawnow;"); // to draw the output to the bash
  
//  mexPrintf("the first triangle is: (%i, %i, %i)\n",t_in[0],t_in[size_t],t_in[2*size_t]);
//  mexPrintf("the second triangle is: (%i, %i, %i)\n",t_in[1],t_in[size_t+1],t_in[2*size_t+1]);
//  mexPrintf("the third triangle is: (%i, %i, %i)\n\n",t_in[2],t_in[size_t+2],t_in[2*size_t+2]);

  for(i=0;i<size_p;i++)
  {
      if (i%1 == 0) 
      {
          mexPrintf("Doing job on point %li of %i\n",i,size_p);
          mexEvalString("drawnow;");
      }
      p_cx = p_in[i];
      p_cy = p_in[size_p+i];
      
      for(iz=0;iz<mesh_z;iz++)
      {
          
          p_cz = iz*z_mesh;
          
          importf(i,iz,importance,N_rays,1);
//		  mexPrintf("Importance Function works");
//          mexEvalString("drawnow;");

          
          for(j=0;j<N_cells;j++)
          {
              // remember that we already corrected the matlab/C indexing issue!
              t_1 = t_in[j];
              t_2 = t_in[N_cells + j];
              t_3 = t_in[2*N_cells + j];
                  
              for(k=0;k<(mesh_z-1);k++)
              {
                              
                  for(l=0;l<N_rays[j+k*N_cells];l++)
                  {
//                        generate the random numbers in the triangle and the z-coordinate
                      u = genrand_real3();
                      v = genrand_real3();
                    
                      if((u+v)>1)
                      {
                          u = 1-u;
                          v = 1-v;
                      }
                  
                       w = 1-u-v;
                      z_rand = (k + genrand_real3())*z_mesh;
                  
                      x_rand = p_in[t_1]*u + p_in[t_2]*v + p_in[t_3]*w;
                      y_rand = p_in[size_p + t_1]*u + p_in[size_p + t_2]*v + p_in[size_p + t_3]*w;                      
//                      random position generation end
                    
//                    here distance is square of distance! see above definition of distance
                      ray_life = ALIVE;
                      
//                      how to include the multiple reflections:
//                      a) generate the mirror point and the corresponding point on the face
//                      b) make transport to the point, then calculate the next point on the next surface
//                      c) multiply all the return gain (propagation), old ones have to be stored in memory

                      gain = beta_v[j+k*N_cells]*propagation(x_rand, y_rand, z_rand, p_cx, p_cy, p_cz, j, k, 1);

                      phi[i+iz*(size_p)] += gain*importance[j+k*N_cells];
//                      return;
                  
                  } // rays loop end
              }//mesh_z loop end
          }//N_cells  loop end
          phi[i+iz*(size_p)] = phi[i+iz*(size_p)]/NumRays;
      }//iz llop end
  }//size_p loop end
  
  mexPrintf("calculations finished, givig back the data\n");
  mexEvalString("drawnow;");
  
  
}// mexfunction end

// write a function for the intergration routine - exp is quite slow
double LineIntegralMCRK4_S(int nSteps, int i_t, int i_z, double PhotonsIn, double length){
	// nSteps	: Number of substeps per cell (more than 5 makes no difference)
    // i_t      : Current triangle
    // i_z      : Current z-level
	// PhotonsIn: Number of photons coming into the first segment to be evaluated
	// length	: Length of the integration path


	double PhotonsOut;
	double h, gain, K1, K2, K3, K4;
	int j;

	h = length / nSteps;
	//printf("%.3e %.3e %.3e %.3e %.3e\n",Photons(ii,jj,kk),VCell,SigmaA,SigmaE,Ntot);
	gain = N_tot*(beta_v[i_t+N_cells*i_z]*(sigma_a+sigma_e)-sigma_a);
	
	//printf("before steps\n");

	// now the RK4 starts
	for (j = 0; j<nSteps; j++){
		// the RK4 by its self
		K1 = PhotonsIn * gain;
		K2 = (PhotonsIn + h/2 * K1)*gain;
		K3 = (PhotonsIn + h/2 * K2)*gain;
		K4 = (PhotonsIn + h * K3)*gain;
		
		//printf("Number of Photons coming in: %.3f - Gain in cell is : %.3f\n",PhotonsIn,gain);

		PhotonsOut = PhotonsIn + h/6*(K1+2*K2+2*K3+K4);

		// RK4 for this step ended
		// now take out the photons and recalculate the gain for saturation
		PhotonsIn=PhotonsOut;
		//printf("PhotonsOut in LineintegralMCRK4_S berechnet zu: %.3E\n",PhotonsOut);
	}// RK4 total end for 1 Cell
	//printf("after steps\n");
//	
	return PhotonsOut;
//    return 0;
//	// be careful with saturation and if the amount of photons substracted is in
//	// bigger than the amount of photons in the cell!!!
}

// write a function for the transport from one point to another 
// as parameters you should give the current position, the destination and the number of reflections wanted
// as well as the beta informations
// give back to integrated gain

double propagation(double x_pos, double y_pos, double z_pos, double x_dest, double y_dest, double z_dest, int t_start, int mesh_start, int N_refl){
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
	int ct; // used to read out cell_type 
    
    
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
		
		
//		  read, which type of cell it is you are propagation in
		ct = cell_type[tri];
        
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
                if (length_help < length && length_help > 0.0)
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
                if (length_help < length && length_help > 0.0)
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
                if (length_help < length && length_help > 0.0)
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
                if (length_help < length && length_help > 0.0)
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
                if (length_help < length && length_help > 0.0)
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

		if (ct == clad_num){
			gain = gain * exp(-clad_abs * length);
		}
		else {
			gain = gain * exp(N_tot*(beta_v[tri+cell_z*N_cells]*(sigma_e + sigma_a)-sigma_a)*length);
		}
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
            ray_life = DEAD;
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

//do the routine which calculates the importance factor and the amount of rays for each of the cells
void importf(int point, int mesh_start, double *importance, int *N_rays, int N_reflections)
{
    int i_t, i_z, Rays_dump=0, rays_left, i_r, rand_t, rand_z;
    double sum_phi=0.0, surf_tot=0.0;
    double distance, x_pos, y_pos, z_pos;
    double prop;
    
//    calculate the gain from the centers of each of the boxes to the observed point
//    calculate the gain and make a "mapping"
//    receipt: pick the point in the center of one cell, 
//    calculate the gain from this point to the observed point,
//    estimate the inner part of the Phi_ASE - Integral,
//    scale the amount of rays proportionally with it
//    sum the amount of rays and scale it to Int=1, which gives the inverse weights
//    the number of rays is determined via floor(), with ceil(), zero-redions could be added
    
//    use the routine "propagation"!, test: no reflections, just exponential
    x_pos = p_in[point];
    y_pos = p_in[point+size_p];
    z_pos = mesh_start * z_mesh;
    
    for (i_t=0;i_t<N_cells;i_t++)
    {
          
        for (i_z=0;i_z<(mesh_z-1);i_z++) //remember the definition differences MatLab/C for indices
        {
//            at this point replace the following routine with propagation(...)
//            later expand this with the beta/tau values...
            prop = propagation(center_x[i_t], center_y[i_t], z_mesh*(i_z+0.5), x_pos, y_pos, z_pos, i_t, i_z, N_reflections);
            importance[i_t + i_z*N_cells] = beta_v[i_t+i_z*N_cells]*(prop);
            sum_phi += importance[i_t + i_z*N_cells];
            
        }
        surf_tot += surface[i_t];
    }
    
//    now calculate the number of rays
    for (i_t=0;i_t<N_cells;i_t++)
    {
        for (i_z=0;i_z<(mesh_z-1);i_z++) //remember the definition differences MatLab/C for indices
        {
//            this is the amount of the sampled rays out of the cells
            N_rays[i_t + i_z*N_cells] = (int)(floor(importance[i_t + i_z*N_cells]/sum_phi*NumRays));
            
            Rays_dump +=  N_rays[i_t + i_z*N_cells];
        }
    }
    
    rays_left = NumRays-Rays_dump;
//    distribute the remaining not distributed rays randomly
    if ((rays_left)>0)
    {
        for (i_r=0;i_r<rays_left;i_r++)
        {
            rand_t = (int )(genrand_real3()*N_cells);
            rand_z = (int )(genrand_real3()*(mesh_z-1));
            N_rays[rand_t + rand_z*N_cells]++;
        }
    }
    
//    now think about the mount of rays which would come out of this volume(surface)
//    dividing this number with the new amount of rays gives the final importance weight for this area!
    for (i_t=0;i_t<N_cells;i_t++)
    {
        for (i_z=0;i_z<(mesh_z-1);i_z++) //remember the definition differences MatLab/C for indices
        {
//            this is the amount of the sampled rays out of the cells
            if (N_rays[i_t + i_z*N_cells]>0)
            {
                importance[i_t + i_z*N_cells] = NumRays*surface[i_t]/surf_tot/N_rays[i_t + i_z*N_cells];
//                importance[i_t + i_z*N_cells] = NumRays*surface[i_t]/surf_tot;
            }
            else
            {
                importance[i_t + i_z*N_cells] = 0; // case of beta of this point == 0 e.g.
            }
        }
    }
}
