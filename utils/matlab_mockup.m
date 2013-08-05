% the folder, where the octrace binary is located
cd ../bin

% the location of the mockup data
load('../utils/testdata_2/save_30.mat')

% the function call to the wrapper script
dndt_ASE = run_octrace(p,normals_x,normals_y,forbidden,normals_p,sorted_int,t_int,z_mesh,mesh_z,N_tot,beta_vol,laser,crystal,beta_cell,surface,x_center,y_center,NumRays);


