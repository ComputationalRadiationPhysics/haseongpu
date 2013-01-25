#ifndef GET_MX_DATA_H
#define GET_MX_DATA_H

host_p_in = (double *)mxGetData(prhs[0]); //point coordinates in 2D . at first size_p x-values, then size_p y-values
host_n_x = (double *)mxGetData(prhs[4]); //normals of the facets, x-components // probably size_t values?? 
host_n_y = (double *)mxGetData(prhs[5]); //normals of the facets, y-components
host_forbidden = (int *)mxGetData(prhs[11]);//are the correspondance of the face used in the previous triangle, which must not be tested (rounding errors)
host_t_in = (int *)mxGetData(prhs[1]);  //association triangle-points - c-indexing-sytle!
host_n_p = (int *)mxGetData(prhs[10]); // gives the Index to one of the points in p_in which is in the plane of the normals - c-indexing-sytle!
host_neighbors = (int *)mxGetData(prhs[6]); //for each cell in t_in, its neighboring cell in plane geometry  - c-indexing-sytle!
host_size_t = (int )mxGetM(prhs[1]); //number of triangles per sheet
N_cells = host_size_t;
size_p = (int )mxGetM(prhs[0]); //number of points
z_mesh = (double)((double *)mxGetData(prhs[14]))[0];
size_z = (int )mxGetN(prhs[2]); //number of meshing in z-direction
host_clad_abs = (double)(cladabs[0]);

#endif
