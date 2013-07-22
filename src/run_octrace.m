%matlab_connector.m

function create_octrace_input (p,normals_x,normals_y,forbidden,normals_p,sorted_int,t_int,z_mesh,mesh_z,N_tot,beta_vol,laser,crystal,beta_cell,surface,x_center,y_center,NumRays)

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

  x=fopen('surface.txt','w')
  fprintf(x,'%.50f\n',surface)
  fclose(x)

  x=fopen('x_center.txt','w')
  fprintf(x,'%.50f\n',x_center)
  fclose(x)

  x=fopen('y_center.txt','w')
  fprintf(x,'%.50f\n',y_center)
  fclose(x)

  x=fopen('NumRays.txt','w')
  fprintf(x,'%d\n',NumRays)
  fclose(x)
end 
  

function [dndt_ASE] = parse_octrace_output ()
  dndt_ASE = load('dndt_ASE.txt');
end

function clean_IO_files ()
  delete('p_in.txt');
  delete('n_x.txt');
  delete('n_y.txt');
  delete('forbidden.txt');
  delete('n_p.txt');
  delete('neighbors.txt');
  delete('t_in.txt');
  delete('z_mesh.txt');
  delete('mesh_z.txt');
  delete('size_t.txt');
  delete('size_p.txt');
  delete('n_tot.txt');
  delete('beta_v.txt');
  delete('sigma_a.txt');
  delete('sigma_e.txt');
  delete('tfluo.txt');
  delete('beta_cell.txt');
  delete('surface.txt');
  delete('x_center.txt');
  delete('y_center.txt');
  delete('dndt_ASE.txt');
  delete('NumRays.txt');
end

function [dndt_ASE] = run_octrace(p,normals_x,normals_y,forbidden,normals_p,sorted_int,t_int,z_mesh,mesh_z,N_tot,beta_vol,laser,crystal,beta_cell,surface,x_center,y_center,NumRays)

  clean_IO_files;

  create_octrace_input(p,normals_x,normals_y,forbidden,normals_p,sorted_int,t_int,z_mesh,mesh_z,N_tot,beta_vol,laser,crystal,beta_cell,surface,x_center,y_center,NumRays);

  %octrace;

  dndt_ASE = parse_octrace_output;
end
