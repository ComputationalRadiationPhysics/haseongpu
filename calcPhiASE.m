% calcPhiASE
% calculates the phi_ASE values for a given input
%first, enter all parameters for the mesh:
% p
% normals_x
% normals_y
% forbidden
% normals_p
% sorted_int
% t_int
% z_mesh 
% mesh_z
% N_tot
% beta_vol
% laser                      struct: parameters for the laser (wavelength sigmas)
% crystal                    struct: contains the crystal parameters (flourescence)
% beta_cell
% surface
% x_center
% y_center
% NumRays                     the minimal number of rays to start for each samplepoint
% clad_int
% clad_number
% clad_abs
% refractiveIndices         the refractive indices for each surface [top_inner, top_outer, bottom_inner, bottom_outer]
% reflectivities            the amount of reflection in case there is no total internal reflection)
% MaxRays                   the maximum number of rays for each samplepoint
% mse_threshold             the maximum allowed expectation (see calc_dndt_ase.cu) for each samplepoint (calculation will be restarted with more rays, if expectation not low enough)
% bool use_reflections      if reflections should be used
%
% phi_ASE                   phi_ASE values, multiple Wavelengths in different columns
% mse_values           real expectation-values for each samplepoint (aligned like phi_ASE)
% N_rays                    number of rays that were used for each samplepoint

%function [phi_ASE, mse_values, N_rays] = calcPhiASE(p,normals_x,normals_y,forbidden,normals_p,sorted_int,t_int,z_mesh,mesh_z,N_tot,beta_vol,laser,crystal,beta_cell,surface,x_center,y_center,NumRays,clad_int, clad_number, clad_abs, refractiveIndices, reflectivities,MaxRays,mse_threshold,use_reflections)

function [phi_ASE, mse_values, N_rays] = calcPhiASE(p,t_int,beta_cell,beta_vol,normals_x,normals_y,sorted_int,surface,x_center,y_center,normals_p,forbidden,NumRays,N_tot,z_mesh,laser,crystal,mesh_z)

%%% Added to the interface %%%
%laser
%crystal
%mesh_z


%%% Defined here in the file %%%
%MaxRays
%mse_threshold
%use_reflections


%%% are not really defined and will be created with dummy values %%% 
%clad_int
%clad_number
%clad_abs
%refractiveIndices
%reflectivities

MaxRays = 100000;
mse=0.05;
use_reflections = false; 
MAX_GPUS=1;
[a,b] = size(p);
minSample=0;
maxSample=(mesh_z*a)-1;
Repetitions=4;

used_dummy = false;

%create dummy variables
if(~exist('mse_threshold','var'))
  %WARNING = 'The variable mse_threshold does not exist'
  [a,b] = size(laser.s_ems);
  mse_threshold = ones(1,a)*mse;
  used_dummy=true;
end

%create dummy variables
if(~exist('clad_int','var') || ~exist('clad_num','var') ||  ~exist('clad_abs','var'))
  %WARNING = 'The variables "clad_int", "clad_num" or "clad_abs" (or a combination) do not exist'
  [a,b] = size(sorted_int);
  clad_int = ones(1,a);
  clad_number = 3;
  clad_abs = 5.5;
  used_dummy=true;
end

%create dummy variables
if(~exist('reflectivities','var') || ~exist('refractiveIndices','var'))
  %WARNING = 'The variables "reflectivities" or "refractiveIndices" (or both) do not exist'
  refractiveIndices = [1.83,1,1.83,1];
  [a,b] = size(sorted_int);
  reflectivities = ones(1,a*2) * 0;
  used_dummy=true;
end

% create the correct reflection-parameter for the c-function
REFLECT='';
if(use_reflections == true)
    REFLECT=' --reflection';
end

if(used_dummy == true)
  disp([ 'WARNING: Some variables were set as dummies' ]);
end

if(~exist('Runmode','var'))
  Runmode='ray_propagation_gpu';
end

CURRENT_DIR = pwd;
FILENAME=[ mfilename('fullpath') '.m' ];
[ CALCPHIASE_DIR, NAME , EXTENSION ] = fileparts(FILENAME);

% create all the textfiles in a separate folder
TMP_FOLDER = [ '/' 'tmp' filesep 'calcPhiASE_tmp' ];

if(strcmpi(Runmode,'mpi'))
  Prefix=[ 'mpiexec -npernode ' num2str(MAX_GPUS) ' ' ];
  % reduce maxGPUs only after setting -npernode
  MAX_GPUS=1;
  % overwrite TMP_FOLDER => needs to be shared among ALL THE NODES!!
  TMP_FOLDER=[ CALCPHIASE_DIR filesep 'mpi_tmp' ];
else
  Prefix='';
end

% make sure that the input is clean 
clean_IO_files(TMP_FOLDER);

% create the new input based on the MATLAB variables
create_calcPhiASE_input(p,normals_x,normals_y,forbidden,normals_p,sorted_int,t_int,z_mesh,mesh_z,N_tot,beta_vol,laser,crystal,beta_cell,surface,x_center,y_center,clad_int,clad_number,clad_abs,refractiveIndices,reflectivities,mse_threshold,TMP_FOLDER,CURRENT_DIR);

  % do the propagation
  system([ Prefix CALCPHIASE_DIR '/bin/calcPhiASE ' '--mode=' Runmode ' --rays=' num2str(NumRays) ' --maxrays=' num2str(MaxRays) REFLECT ' --input=' TMP_FOLDER ' --output=' TMP_FOLDER ' --min_sample_i=' num2str(minSample) ' --max_sample_i=' num2str(maxSample) ' --maxgpus=' num2str(MAX_GPUS) ' --repetitions=' num2str(Repetitions) ]);

  % get the result
  [ mse_values, N_rays, phi_ASE ] = parse_calcPhiASE_output(TMP_FOLDER,CURRENT_DIR);

  % cleanup
  %clean_IO_files(TMP_FOLDER);
  % 
end

%takes all the variables and puts them into textfiles, so the CUDA code can parse them
function create_calcPhiASE_input (p,normals_x,normals_y,forbidden,normals_p,sorted_int,t_int,z_mesh,mesh_z,N_tot,beta_vol,laser,crystal,beta_cell,surface,x_center,y_center,clad_int,clad_number,clad_abs,refractiveIndices,reflectivities,mse_threshold,FOLDER,CURRENT_DIR)

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
fprintf(x,'%d\n',clad_number);
fclose(x);

x=fopen('clad_abs.txt','w');
fprintf(x,'%.50f\n',clad_abs);
fclose(x);

x=fopen('refractive_indices.txt','w');
fprintf(x,'%3.5f\n',refractiveIndices);
fclose(x);

x=fopen('reflectivities.txt','w');
fprintf(x,'%.50f\n',reflectivities);
fclose(x);

x=fopen('mse_threshold.txt','w');
fprintf(x,'%.15f\n',mse_threshold);
fclose(x);

cd(CURRENT_DIR);
end 


% takes the output from the CUDA code and fills it into a variable
function [mseValues,  N_rays, phi_ASE] = parse_calcPhiASE_output (FOLDER,CURRENT_DIR)
cd (FOLDER);
fid = fopen('phi_ASE.txt');
arraySize = str2num(fgetl(fid));
phi_ASE = str2num(fgetl(fid));
phi_ASE = reshape(phi_ASE,arraySize);
fclose(fid);

fid=fopen('mse_values.txt');
arraySize=str2num(fgetl(fid));
mseValues = str2num(fgetl(fid));
mseValues = reshape(mseValues,arraySize);
fclose(fid);


fid = fopen('N_rays.txt');
arraySize = str2num(fgetl(fid));
N_rays = str2num(fgetl(fid));
N_rays = reshape(N_rays,arraySize);
fclose(fid);

cd (CURRENT_DIR);
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

  warning(s);
end

