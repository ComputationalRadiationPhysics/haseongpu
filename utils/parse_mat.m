function parse_mat (folder,file)
format long
save_precision(50)
cd(folder)
load(file)

x=fopen("p_in.txt","w")
fprintf(x,'%.50f\n',p)
fclose(x)

x=fopen("n_x.txt","w")
fprintf(x,'%.50f\n',normals_x)
fclose(x)

x=fopen("n_y.txt","w")
fprintf(x,'%.50f\n',normals_y)
fclose(x)

x=fopen("forbidden.txt","w")
fprintf(x,'%d\n',forbidden)
fclose(x)

x=fopen("n_p.txt","w")
fprintf(x,'%d\n',normals_p)
fclose(x)

x=fopen("neighbors.txt","w")
fprintf(x,'%d\n',sorted_int)
fclose(x)

x=fopen("t_in.txt","w")
fprintf(x,'%d\n',t_int)
fclose(x)

% thickness of one slice!
x=fopen("z_mesh.txt","w")
fprintf(x,'%.50f\n',z_mesh)
fclose(x)

% number of slices
x=fopen("mesh_z.txt","w")
fprintf(x,'%d\n',mesh_z)
fclose(x)


x=fopen("size_t.txt","w")
fprintf(x,'%d\n',rows(t_int))
fclose(x)

x=fopen("size_p.txt","w")
fprintf(x,'%d\n',rows(p))
fclose(x)

% fixed
x=fopen("clad_abs.txt","w")
fprintf(x,'5.5\n')
fclose(x)

%fixed
x=fopen("clad_num.txt","w")
fprintf(x,'3\n')
fclose(x)

x=fopen("n_tot.txt","w")
fprintf(x,'%.50f\n',N_tot)
fclose(x)

x=fopen("beta_v.txt","w")
fprintf(x,'%.50f\n',beta_vol)
fclose(x)

x=fopen("sigma_a.txt","w")
fprintf(x,'%.50f\n',laser.s_abs)
fclose(x)

x=fopen("sigma_e.txt","w")
fprintf(x,'%.50f\n',laser.s_ems)
fclose(x)

x=fopen("tfluo.txt","w")
fprintf(x,'%.50f\n',crystal.tfluo)
fclose(x)

% fixed
x=fopen("cell_type.txt","w")
fprintf(x,'%d\n',ones(1,rows(sorted_int)))
fclose(x)


%x=fopen("NumRays.txt","w")
%fprintf(x,'%.50f\n',NumRays)
%fclose(x)

x=fopen("beta_cell.txt","w")
fprintf(x,'%.50f\n',beta_cell)
fclose(x)

%x=fopen("surface.txt","w")
%fprintf(x,'%.50f\n',surface)
%fclose(x)

%x=fopen("x_center.txt","w")
%fprintf(x,'%.50f\n',x_center)
%fclose(x)

%x=fopen("y_center.txt","w")
%fprintf(x,'%.50f\n',y_center)
%fclose(x)

end 
