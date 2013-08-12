% run_octrace
% calculates the dndt_ASE values for a given input
function [expected_values, importance, N_rays, phi_ASE] = run_octrace(p,normals_x,normals_y,forbidden,normals_p,sorted_int,t_int,z_mesh,mesh_z,N_tot,beta_vol,laser,crystal,beta_cell,surface,x_center,y_center,NumRays,clad_int, clad_number, clad_abs)

  % create all the textfiles in a separate folder
  TMP_FOLDER = 'octrace_tmp';
  FOLDER = [ pwd filesep TMP_FOLDER ];

  % make sure that the input is clean 
  clean_IO_files(FOLDER);

  % create the new input based on the MATLAB variables
  create_octrace_input(p,normals_x,normals_y,forbidden,normals_p,sorted_int,t_int,z_mesh,mesh_z,N_tot,beta_vol,laser,crystal,beta_cell,surface,x_center,y_center,clad_int,clad_number,clad_abs,FOLDER);

  % do the propagation
  system(['./octrace ' '--mode=ray_propagation_gpu ' '--silent ' '--rays=' num2str(NumRays) ' --experiment=' FOLDER ]);

  % get the result
  [ expected_values, importance, N_rays, phi_ASE ] = parse_octrace_output;

  % cleanup
  clean_IO_files(FOLDER);
end

%takes all the variables and puts them into textfiles, so the CUDA code can parse them
function create_octrace_input (p,normals_x,normals_y,forbidden,normals_p,sorted_int,t_int,z_mesh,mesh_z,N_tot,beta_vol,laser,crystal,beta_cell,surface,x_center,y_center,clad_int,clad_number,clad_abs,FOLDER)

  mkdir(FOLDER);
  cd(FOLDER);

  x=fopen('p_in.txt','w');
  fprintf(x,'%.50f\n',p);
  fclose(x);

  x=fopen('n_x.txt','w');
  fprintf(x,'%.50f\n',normals_x);
  fclose(x);

  x=fopen('n_y.txt','w');
  fprintf(x,'%.50f\n',normals_y);
  fclose(x);

  x=fopen('forbidden.txt','w');
  fprintf(x,'%d\n',forbidden);
  fclose(x);

  x=fopen('n_p.txt','w');
  fprintf(x,'%d\n',normals_p);
  fclose(x);

  x=fopen('neighbors.txt','w');
  fprintf(x,'%d\n',sorted_int);
  fclose(x);

  x=fopen('t_in.txt','w');
  fprintf(x,'%d\n',t_int);
  fclose(x);

  % thickness of one slice!
  x=fopen('z_mesh.txt','w');
  fprintf(x,'%.50f\n',z_mesh);
  fclose(x);

  % number of slices
  x=fopen('mesh_z.txt','w');
  fprintf(x,'%d\n',mesh_z);
  fclose(x);

  x=fopen('size_t.txt','w');
  [a,b] = size(t_int);
  fprintf(x,'%d\n',a);
  fclose(x);

  x=fopen('size_p.txt','w');
  [a,b] = size(p);
  fprintf(x,'%d\n',a);
  fclose(x);

  x=fopen('n_tot.txt','w');
  fprintf(x,'%.50f\n',N_tot);
  fclose(x);

  x=fopen('beta_v.txt','w');
  fprintf(x,'%.50f\n',beta_vol);
  fclose(x);

  x=fopen('sigma_a.txt','w');
  fprintf(x,'%.50f\n',laser.s_abs);
  fclose(x);

  x=fopen('sigma_e.txt','w');
  fprintf(x,'%.50f\n',laser.s_ems);
  fclose(x);

  x=fopen('tfluo.txt','w');
  fprintf(x,'%.50f\n',crystal.tfluo);
  fclose(x);

  x=fopen('beta_cell.txt','w');
  fprintf(x,'%.50f\n',beta_cell);
  fclose(x);

  x=fopen('surface.txt','w');
  fprintf(x,'%.50f\n',surface);
  fclose(x);

  x=fopen('x_center.txt','w');
  fprintf(x,'%.50f\n',x_center);
  fclose(x);

  x=fopen('y_center.txt','w');
  fprintf(x,'%.50f\n',y_center);
  fclose(x);

  x=fopen('clad_int.txt','w');
  fprintf(x,'%d\n',clad_int);
  fclose(x);

  x=fopen('clad_num.txt','w');
  fprintf(x,'%d\n',clad_num);
  fclose(x);

  x=fopen('clad_abs.txt','w');
  fprintf(x,'%.50f\n',clad_abs);
  fclose(x);

  cd ..
end 
  

% takes the output from the CUDA code and fills it into a variable
function [expectedValues, importance, N_rays, phi_ASE] = parse_octrace_output ()
  phi_ASE = load('phi_ASE.txt');
  expectedValues = load('expected_values.txt');
  importance = load('importance.txt');
  N_rays = load('N_rays.txt');
end

% deletes the temporary folder and the dndt_ASE.txt
function clean_IO_files (TMP_FOLDER)
  
  % disable warnings for this (if file is nonexistant...)
  s = warning;
  warning off all;

  A = exist(TMP_FOLDER,'dir');

  if A == 7
    rmdir(TMP_FOLDER,'s');
  end

  delete('dndt_ASE.txt');
  delete('expected_values.txt');
  delete('importance.txt');
  delete('N_rays.txt');

  warning(s);
end

