
% clear all

% time testing the for_loops routine
% load variable.mat

% run ASE_3D_triangles and stopp it before the computation starts, this
% give all the necessary variables!

NumRays = 1e5;
NumRays = int32(NumRays);

tic
[rand_array, phi_ASE, importance, N_rays] = for_loops(p,t_int,beta_cell,beta_vol,normals_x,normals_y,sorted_int,surface,x_center,y_center,normals_p,forbidden, NumRays);
toc

surface_total = sum(surfaces);
volume_total = surface_total*D;

phi_ASE = phi_ASE./(4*3.1415)*volume_total;

% now plot it and see if the random variables lie withing the rectangle or
% not

%% use last triangle
% % Px(1) = p(t_int(N_cells,1),1);
% % Py(1) = p(t_int(N_cells,1),2);
% % 
% % Px(2) = p(t_int(N_cells,2),1);
% % Py(2) = p(t_int(N_cells,2),2);
% % 
% % Px(3) = p(t_int(N_cells,3),1);
% % Py(3) = p(t_int(N_cells,3),2);
% N_test = 1;
% Px(1) = p(t_int(N_test,1)+1,1);
% Py(1) = p(t_int(N_test,1)+1,2);
% 
% Px(2) = p(t_int(N_test,2)+1,1);
% Py(2) = p(t_int(N_test,2)+1,2);
% 
% Px(3) = p(t_int(N_test,3)+1,1);
% Py(3) = p(t_int(N_test,3)+1,2);
% 
% X = [Px,Px(1)];
% Y = [Py,Py(1)];
% 
% rand_x = 0.737855;
% rand_y = 1.302254;
% 
% length = 0.437257;%0.062584;
% 
% vec_x = -0.796131;
% vec_y = -0.605124;
% 
% pos_x = rand_x + length*vec_x;
% pos_y = rand_y + length*vec_y;
% 
% % now show the triangle and returned random positions
% hold on
% plot(X,Y);
% % plot(rand_array(:,1),rand_array(:,2),'*');
% plot(rand_x,rand_y,'*');
% plot(pos_x,pos_y,'*');
% hold off
% % 
% % %% final plot
% % trimesh(double(t_int),double(p(:,1)),double(p(:,2)),phi_ASE(:,1))
% % view(2),axis equal,axis off,drawnow

%% new plot with interpolation
figure(2)
[x_grid,y_grid]=meshgrid(-1.5:0.01:1.5);
Phi_ASE_Interp = griddata(p(:,1),p(:,2),phi_ASE(:,1),x_grid,y_grid);
imagesc(Phi_ASE_Interp);
axis equal;
colorbar;