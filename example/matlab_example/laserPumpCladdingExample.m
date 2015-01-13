 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Copyright 2013 Daniel Albach, Erik Zenker, Carlchristian Eckert
 %
 % This file is part of HASEonGPU
 %
 % HASEonGPU is free software: you can redistribute it and/or modify
 % it under the terms of the GNU General Public License as published by
 % the Free Software Foundation, either version 3 of the License, or
 % (at your option) any later version.
 %
 % HASEonGPU is distributed in the hope that it will be useful,
 % but WITHOUT ANY WARRANTY; without even the implied warranty of
 % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 % GNU General Public License for more details.
 %
 % You should have received a copy of the GNU General Public License
 % along with HASEonGPU.
 % If not, see <http://www.gnu.org/licenses/>.
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Laser pump routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASE calulations main file 16kW
% longer calculation without pump action possibel for heat generation
% estimation in the cladding
%
% @author Daniel Albach
% @date   2011-05-03
% @licence: GPLv3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
tic
isOctave = exist('OCTAVE_VERSION') ~= 0;
if (isOctave)
   page_output_immediately(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Definitions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Crystal parameter 
crystal.doping = 2;
crystal.length = 0.7;    % [cm]
crystal.tfluo = 9.5e-4;  % 1.0e-3 pour DA
crystal.nlexp = 1;
crystal.levels = 10;

% Timesteps
steps.time = 100;
steps.crys = crystal.levels;

% Pump parameter
pump_stretch = 1;
pump_aspect = 1;
pump_diameter = 3;
pump.s_abs = 0.778e-20;                % Absorption cross section in cm² (0.778e-20 pour DA)
pump.s_ems = 0.195e-20;                % Emission cross section in cm²(0.195e-20 pour DA)
pump.I = 16e3;                         % Pump intensity in W/cm²
pump.T = 1e-3;                         % Pump duration in s
pump.wavelength = 940e-9;              % Pump wavelength in m
pump.ry = pump_diameter/2*pump_aspect;
pump.rx = pump_diameter/2;             % Pump radius
pump.exp = 40;                         % Pump exponent

% Laser parameter
laser.s_abs = load('sigma_a.txt');      % Absorption spectrum cm2(1.16e-21 pour DA)
laser.s_ems = load('sigma_e.txt');      % Emission spectrum in cm2(2.48e-20)
laser.l_abs = load('lambda_a.txt');     % Wavelengths absoroption spectrum in nm (x values for absorption)
laser.l_ems = load('lambda_e.txt');     % Wavelengths emission spectrum in nm (y values for emission)
laser.l_res = 1000;                     % Resolution of linear interpolated spectrum
laser.I = 1e6;                          % Laser intensity
laser.T = 1e-8;                         % Laser duration
laser.wavelength = 1030e-9;             % Lasing wavelength in m
[laser.max_ems, i] = max(laser.s_ems);
laser.max_abs = laser.s_abs(i);

% Mode parameter
mode.BRM = 1;  % 1 Backreflection, 0 Transmissionmode
mode.R = 1;    % reflectivity of the coating
mode.extr =0;  %no extraction!!

% Constants
const.N1per = 1.38e20;
const.c = 3e8;
const.h = 6.626e-34;
N_tot = const.N1per*crystal.doping;
Ntot_gradient = zeros(crystal.levels, 1);      % doping gradient if not given by load!
Ntot_gradient(:) = crystal.doping*const.N1per;
mesh_z = crystal.levels;
z_mesh = crystal.length/(mesh_z-1);
timeslice = 50;
timeslice_tot = 150;
timetotal = 1e-3; %[s]
time_t = timetotal/timeslice;

% ASE application
maxGPUs           = 2;
nPerNode          = 4;
deviceMode             = 'gpu';
%parallelMode           = 'mpi';
parallelMode          = 'threaded';
useReflections    = true;
refractiveIndices = [1.83,1,1.83,1];
repetitions       = 4;
minRaysPerSample  = 1e5;
maxRaysPerSample  = minRaysPerSample * 100;
mseThreshold      = 0.005;


% Constants for short use
c = const.c; %m/s
h = const.h; %Js

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the grid from file
% the grid is done using the distmesh routines
% use e.g. cr_60mm_30mm.m in meshing folder
load('pt.mat');

% Create points and triangles from grid
set_variables(p,t);
load('variable.mat');

% Mesh dependand definitions
N_cells = size(t,1);
beta_cell = zeros(size(p,1),mesh_z);
beta_vol = zeros(N_cells,mesh_z-1);
[nT,b]            = size(sorted_int);
reflectivities    = zeros(1,nT*2);

% In first timeslice phi_ASE stays 0
phi_ASE=zeros(size(p,1),mesh_z);
dndt_ASE(:,:) = phi_ASE(:,:);
flux_clad = phi_ASE(:,:);
dndt_pump = (beta_cell);

% Save values in file -> first timestep the phi_ASE is zero
vtk_wedge('beta_cell_0.vtk',beta_cell, p, t_int, mesh_z, z_mesh);
vtk_wedge('dndt_pump_0.vtk',dndt_pump, p, t_int, mesh_z, z_mesh);
vtk_wedge('dndt_ASE_0.vtk',dndt_ASE , p, t_int, mesh_z, z_mesh);
save('save_0.mat');

temp = pump.T;
time_beta = 1e-6;
pump.T = time_beta;
temp_f = crystal.tfluo;
crystal.tfluo = 1;
beta_c_2 = zeros(size(p,1),mesh_z);
intensity = pump.I;

% Now call beta_int for each of the points
for i_p=1:size(p,1)
   beta_crystal=beta_cell(i_p,:);
   pulse = zeros(steps.time,1);
   pump.I = intensity*exp(-sqrt(p(i_p,1)^2/pump.ry^2+p(i_p,2)^2/pump.rx^2)^pump.exp);
   [beta_crystal,beta_store,pulse,Ntot_gradient] = beta_int3(beta_crystal,pulse,const,crystal,steps,pump,mode,Ntot_gradient);
   beta_c_2(i_p,:)=beta_crystal(:);
end

% don't foget, that we pump in -z-direction, but the geometry direction id
% inverse! So flip the geometry upside-down or define the surface
% conditions different!

dndt_pump = (beta_c_2 - beta_cell)./time_beta;

pump.I = intensity;
pump.T = temp;

crystal.tfluo = temp_f;

time_int = 1e-6;

% Integrate for the first time
for i_p=1:size(p,1)
    for i_z=1:mesh_z
        beta_cell(i_p,i_z) = crystal.tfluo*(dndt_pump(i_p,i_z)-dndt_ASE(i_p,i_z))*(1-exp(-time_t/crystal.tfluo)) + beta_cell(i_p,i_z)*exp(-time_t/crystal.tfluo);
    end
end

vtk_wedge('beta_cell_1.vtk',beta_cell, p, t_int, mesh_z, z_mesh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Definition of cladding %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Which index is the cladding?
clad_number           = 1;
clad_number           = int32(clad_number);
[numberOfTriangles,b] = size(t_int);
clad                  = zeros(numberOfTriangles,1);

% Absorption of the cladding
clad_abs = 5.5;

% Prepare which part is concerning the cladding and which one is the Yb doped part
% this part mus be adapted for each case!
clad_index = find(clad==clad_number);
Yb_index = find(clad~=clad_number);

% the logic is inverse, which index to delete ;)
% this represents the triangle index and using the p-variable, you get the points...
Yb_points_t = t;
Yb_points_t(clad_index,:) = [];
Yb_points_t2 = reshape(Yb_points_t,[],1);
clear Yb_points_t;
Yb_points_t = unique(Yb_points_t2);
clear Yb_points_t2;

clad_points_t = t;
clad_points_t(Yb_index,:) = [];
clad_points_t2 = reshape(clad_points_t,[],1);
clear clad_points_t;
clad_points_t = unique(clad_points_t2);
clear clad_points_t2;

% Now get the lengths of the 
length_clad = size(clad_points_t,1);
length_Yb = size(Yb_points_t,1);

% Make sure that clad is integer variable
clad_int = int32(clad);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main pump loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_slice=1:timeslice_tot-1
    disp(['TimeSlice ' num2str(i_slice) ' of ' num2str(timeslice_tot-1) ' started']);
    % ******************* BETA PUMP TEST ******************************
    % make a test with the gain routine "gain.m" for each of the points
    % define the intensity at the nominal points of the grid and make the
    % pumping to get the temporal evolution of the beta
    % make a rectangular field, make the interpolation of the points on it
    % if we don't have any formula discribing the intensity distribution
    % as first test do a supergaussian distribution for 1µs and then estimate
    % the dndt|pump
    
    % Write output for this timeslice to file
    file_b = ['beta_cell_' num2str(i_slice) '.vtk'];
    file_p = ['dndt_pump_' num2str(i_slice) '.vtk'];
    file_A = ['dndt_ASE_' num2str(i_slice) '.vtk'];
    file_C = ['flux_clad_' num2str(i_slice) '.vtk'];
    
    vtk_wedge(file_b,beta_cell, p, t_int, mesh_z, z_mesh);
    vtk_wedge(file_p,dndt_pump, p, t_int, mesh_z, z_mesh);
    
    % Interpolate beta_vol from beta_cell
    x_1 = p(:,1);
    y_1 = p(:,2);
    x_2 = p(:,1);
    y_2 = p(:,2);

    x = vertcat(x_1,x_2);
    y = vertcat(y_1,y_2);

    z_1 = zeros(size(x_1,1),1);
    z_2 = zeros(size(x_2,1),1);
    z_2 = z_2+z_mesh;

    z = vertcat(z_1,z_2);

    xi = x_center;
    yi = y_center;
    zi = zeros(size(xi,1),1);
    zi = zi + z_mesh/2;

    for i_z=1:(mesh_z-1)
        v_1 = beta_cell(:,i_z);
        v_2 = beta_cell(:,i_z+1);

        v = vertcat(v_1,v_2);

        if (isOctave)
          beta_vol(:,i_z) = griddata3(x,y,z,v,xi,yi,zi);
        else
          beta_vol(:,i_z) = griddata(x,y,z,v,xi,yi,zi);
        end

        z = z + z_mesh;
        zi = zi + z_mesh;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Call external ASE application %%%%%%%%%%%%%%%%%%%%
     [phi_ASE, mse_values, N_rays] = calcPhiASE(...
       p,...
       t_int,...
       beta_cell,...
       beta_vol,...
       clad_int,...
       clad_number,...
       clad_abs,...
       useReflections,...
       refractiveIndices,...
       reflectivities,...
       normals_x,...
       normals_y,...
       sorted_int,...
       surface,...
       x_center,...
       y_center,...
       normals_p,...
       forbidden,...
       minRaysPerSample,...
       maxRaysPerSample,...
       mseThreshold,...
       repetitions,...
       N_tot,...
       z_mesh,...
       laser,...
       crystal,...
       mesh_z,...
       deviceMode,...
       parallelMode,...
       maxGPUs,...
       nPerNode);


    % Calc dn/dt ASE, change it to only the Yb:YAG points
    for i_p=1:length_Yb
        for i_z=1:mesh_z
            pos_Yb = Yb_points_t(i_p);
	    % Calc local gain (g_l)
            g_l = -(N_tot*laser.max_abs - N_tot*beta_cell(pos_Yb,i_z)*(laser.max_ems+laser.max_abs));
            dndt_ASE(pos_Yb,i_z) = g_l * phi_ASE(pos_Yb,i_z) / crystal.tfluo;
        end
    end
    

    
    % Now for the cladding points
    for i_p=1:length_clad
        for i_z=1:mesh_z
            % Calc local gain (g_l)
            pos_clad = clad_points_t(i_p);
            flux_clad(pos_clad,i_z) = clad_abs*phi_ASE(pos_clad,i_z)/crystal.tfluo;
        end
    end
    
    vtk_wedge(file_A,dndt_ASE , p, t_int, mesh_z, z_mesh);
    vtk_wedge(file_C,flux_clad , p, t_int, mesh_z, z_mesh);
    save(['save_' num2str(i_slice) '.mat']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Prepare next timeslice %%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp = pump.T;
    time_beta = 1e-6;
    pump.T = time_beta;
    temp_f = crystal.tfluo;
    crystal.tfluo = 1;
    beta_c_2 = zeros(size(p,1),mesh_z);
    intensity = pump.I;

    % Now call beta_int for each of the points
    for i_p=1:size(p,1)
        beta_crystal=beta_cell(i_p,:);
        pulse = zeros(steps.time,1);
        pump.I = intensity*exp(-sqrt(p(i_p,1)^2/pump.ry^2+p(i_p,2)^2/pump.rx^2)^pump.exp);
        [beta_crystal,beta_store,pulse,Ntot_gradient] = beta_int3(beta_crystal,pulse,const,crystal,steps,pump,mode,Ntot_gradient);
        beta_c_2(i_p,:) = beta_crystal(:);
    end

    if (i_slice<timeslice)
        dndt_pump = (beta_c_2 - beta_cell)./time_beta;
    else
        dndt_pump = zeros(size(p,1),mesh_z);
    end
    
    pump.I = intensity;
    pump.T = temp;
    crystal.tfluo = temp_f;

    % Recalculate beta_cell
    for i_p=1:size(p,1)
        for i_z=1:mesh_z
            beta_cell(i_p,i_z) = crystal.tfluo*(dndt_pump(i_p,i_z)-dndt_ASE(i_p,i_z))*(1-exp(-time_t/crystal.tfluo))+beta_cell(i_p,i_z)*exp(-time_t/crystal.tfluo);
        end
    end
end

file_b = ['beta_cell_' num2str(timeslice_tot) '.vtk'];
file_p = ['dndt_pump_' num2str(timeslice_tot) '.vtk'];
file_A = ['dndt_ASE_' num2str(timeslice_tot) '.vtk'];
file_C = ['flux_clad_' num2str(timeslice_tot) '.vtk'];

vtk_wedge(file_b,beta_cell, p, t_int, mesh_z, z_mesh);
vtk_wedge(file_p,dndt_pump, p, t_int, mesh_z, z_mesh);
vtk_wedge(file_A,dndt_ASE , p, t_int, mesh_z, z_mesh);
vtk_wedge(file_C,flux_clad , p, t_int, mesh_z, z_mesh);

disp('Calculations finished');
toc

% Obviously I forgot to save the latest results, must be added!
