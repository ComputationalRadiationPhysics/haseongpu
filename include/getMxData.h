#ifndef GET_MX_DATA_H
#define GET_MX_DATA_H

#include <matrix.h>
#include <mex.h>

double p_cx,p_cy, p_cz;

//  mex inputs and outputs
//  double *p_in, *beta_p, *beta_v, *n_x, *n_y, *surface, *center_x, *center_y;
//  int *t_in, *neighbors, *n_p;

mxArray *rand_out, *phi_ASE, *import, *N_r;
double *rand, *phi, *importance;
int size_t;
int *N_rays, *Rays;
double *dop, *zm;
double *sa, *se; // doping information input
double *cladabs; // cladding information
int *cladnum; //cladding information
int N_refl; // number of reflections that the beam will do

double *p_in, *beta_p, *beta_v, *n_x, *n_y, *surface, *center_x, *center_y;
int *t_in, *n_p, *neighbors, *forbidden;
int size_p, N_cells, mesh_z, rays;
int NumRays;
double z_mesh=0.0356;
double clad_abs;
int clad_num, *cell_type;

//  spectral properties
double beta = 0.12;
double sigma_a = 1.16e-21;
double sigma_e = 2.48e-20;
double doping = 5;
double N_tot = 1.38e20;


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

#endif
