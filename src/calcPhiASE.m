%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calcPhiASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% calculates the phiASE values for a given input
% most meshing paramers are given through the function parameters.
% However, many parameters for optimization of the computation are
% set in the beginning of the function (adjust as needed)
% 
% for most mesh parameters see README file
%
% @return phiASE the ASE-Flux for all the given sample points
% @return mseValues the MeanSquaredError values corresponding to phiASE
% @return raysUsedPerSample the number of rays used to calculate each phiASE value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phiASE, mseValues, raysUsedPerSample] = calcPhiASE(points,trianglePointIndices,betaCells,betaVolume,triangleNormalsX,triangleNormalsY,triangleNeighbors,triangleSurfaces,triangleCenterX,triangleCenterY,triangleNormalPoint,forbiddenEdge,minRaysPerSample,nTot,thickness,laserParameter,crystal,numberOfLevels)

%%%%% configuration of values that are not reachable by interface %%%%%%%%%
maxGPUs=1;  % should be 1 for use with MPI, otherwise the number of devices in the node
nPerNode=4; % for a MPI node with 4 devices
Runmode='mpi';
%Runmode='ray_propagation_gpu'; % runmode for the threaded approach

maxRaysPerSample = 100000000;
mseThreshold=0.05;
Repetitions=4;

useReflections = true; 
refractiveIndices = [1.83,1,1.83,1];
[nT,b] = size(triangleNeighbors);
reflectivities = zeros(1,nT*2);
claddingInt = ones(1,nT);
claddingNumber = 3;
claddingAbsorption = 5.5;


%%%%%%%% overwriting values from the interface %%%%%%%%%%%%
minRaysPerSample = 10000000;


%%%%%%%%%%%%% auto-generating some more input %%%%%%%%%%%%%
minSample=0;
[nP,b] = size(points);
maxSample=(numberOfLevels*nP)-1;

if(useReflections == true)
    REFLECT=' --reflection';
else
    REFLECT='';
end

if(strcmpi(Runmode,'mpi'))
  Prefix=[ 'mpiexec -npernode ' num2str(nPerNode) ' ' ];
  maxGPUs=1;
else
  Prefix='';
end

% create a tmp-folder in the same directory as this script
FILENAME=[ mfilename('fullpath') '.m' ];
[ CALCPHIASE_DIR, NAME , EXTENSION ] = fileparts(FILENAME);
TMP_FOLDER=[ CALCPHIASE_DIR filesep 'input_tmp' ];


%%%%%%%%%%%%%%%%%% doing the computation %%%%%%%%%%%%%%%%%%
% make sure that the temporary folder is clean 
clean_IO_files(TMP_FOLDER);

% create new input in the temporary folder
create_calcPhiASE_input(points,triangleNormalsX,triangleNormalsY,forbiddenEdge,triangleNormalPoint,triangleNeighbors,trianglePointIndices,thickness,numberOfLevels,nTot,betaVolume,laserParameter,crystal,betaCells,triangleSurfaces,x_center,triangleCenterY,claddingInt,claddingNumber,claddingAbsorption,refractiveIndices,reflectivities,TMP_FOLDER);

% do the propagation
status = system([ Prefix CALCPHIASE_DIR '/bin/calcPhiASE ' '--mode=' Runmode ' --rays=' num2str(minRaysPerSample) ' --maxrays=' num2str(maxRaysPerSample) REFLECT ' --input=' TMP_FOLDER ' --output=' TMP_FOLDER ' --min_sample_i=' num2str(minSample) ' --max_sample_i=' num2str(maxSample) ' --maxgpus=' num2str(maxGPUs) ' --repetitions=' num2str(Repetitions) ' --mse-threshold=' num2str(mseThreshold) ]);

if(status ~= 0)
    error(['this step of the raytracing computation did NOT finish successfully. Aborting.']);
end

% get the result
[ mseValues, raysUsedPerSample, phiASE ] = parse_calcPhiASE_output(TMP_FOLDER);

% cleanup
clean_IO_files(TMP_FOLDER);

end





%%%%%%%%%%%%%%%%%%%%%%%%% parse_calcPhiASE_output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% takes all the variables and puts them into textfiles
% so the CUDA code can parse them. Take care that the
% names of the textfiles match those in the parsing function
% 
% for most parameters, see calcPhiASE (above)
% @param FOLDER the folder in which to create the input files (usually a
%               temporary folder visible by all the participating nodes)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function create_calcPhiASE_input (points,triangleNormalsX,triangleNormalsY,forbiddenEdge,triangleNormalPoint,triangleNeighbors,trianglePointIndices,thickness,numberOfLevels,nTot,betaVolume,laserParameter,crystal,betaCells,triangleSurfaces,triangleCenterX,triangleCenterY,claddingInt,claddingNumber,claddingAbsorption,refractiveIndices,reflectivities,FOLDER)
CURRENT_DIR = pwd;

mkdir(FOLDER);
cd(FOLDER);

%x=fopen('points.txt','w');
x=fopen('p_in.txt','w');
fprintf(x,'%.50f\n',points);
fclose(x);

%x=fopen('triangleNormalsX.txt','w');
x=fopen('n_x.txt','w');
fprintf(x,'%.50f\n',triangleNormalsX);
fclose(x);

%x=fopen('triangleNormalsY.txt','w');
x=fopen('n_y.txt','w');
fprintf(x,'%.50f\n',triangleNormalsY);
fclose(x);

%x=fopen('forbiddenEdge.txt','w');
x=fopen('forbidden.txt','w');
fprintf(x,'%d\n',forbiddenEdge);
fclose(x);

%x=fopen('triangleNormalPoint.txt','w');
x=fopen('n_p.txt','w');
fprintf(x,'%d\n',triangleNormalPoint);
fclose(x);

%x=fopen('triangleNeighbors.txt','w');
x=fopen('neighbors.txt','w');
fprintf(x,'%d\n',triangleNeighbors);
fclose(x);

%x=fopen('trianglePointIndices.txt','w');
x=fopen('t_in.txt','w');
fprintf(x,'%d\n',trianglePointIndices);
fclose(x);

% thickness of one slice!
%x=fopen('thickness.txt','w');
x=fopen('z_mesh.txt','w');
fprintf(x,'%.50f\n',thickness);
fclose(x);

% number of slices
%x=fopen('numberOfLevels.txt','w');
x=fopen('mesh_z.txt','w');
fprintf(x,'%d\n',numberOfLevels);
fclose(x);

%x=fopen('numberOfTriangles.txt','w');
x=fopen('size_t.txt','w');
[a,b] = size(trianglePointIndices);
fprintf(x,'%d\n',a);
fclose(x);

%x=fopen('numberOfPoints','w');
x=fopen('size_p.txt','w');
[a,b] = size(points);
fprintf(x,'%d\n',a);
fclose(x);

%x=fopen('nTot.txt','w');
x=fopen('n_tot.txt','w');
fprintf(x,'%.50f\n',nTot);
fclose(x);

%x=fopen('betaVolume.txt','w');
x=fopen('beta_v.txt','w');
fprintf(x,'%.50f\n',betaVolume);
fclose(x);

%x=fopen('sigmaA.txt','w');
x=fopen('sigma_a.txt','w');
fprintf(x,'%.50f\n',laserParameter.s_abs);
fclose(x);

%x=fopen('sigmaE.txt','w');
x=fopen('sigma_e.txt','w');
fprintf(x,'%.50f\n',laserParameter.s_ems);
fclose(x);

%x=fopen('tFluo.txt','w');
x=fopen('tfluo.txt','w');
fprintf(x,'%.50f\n',crystal.tfluo);
fclose(x);

%x=fopen('betaCells.txt','w');
x=fopen('beta_cell.txt','w');
fprintf(x,'%.50f\n',betaCells);
fclose(x);

%x=fopen('triangleSurfaces.txt','w');
x=fopen('surface.txt','w');
fprintf(x,'%.50f\n',triangleSurfaces);
fclose(x);

%x=fopen('triangleCenterX.txt','w');
x=fopen('x_center.txt','w');
fprintf(x,'%.50f\n',triangleCenterX);
fclose(x);

%x=fopen('triangleCenterY.txt','w');
x=fopen('y_center.txt','w');
fprintf(x,'%.50f\n',triangleCenterY);
fclose(x);

%x=fopen('claddingInt.txt','w');
x=fopen('clad_int.txt','w');
fprintf(x,'%d\n',claddingInt);
fclose(x);

%x=fopen('claddingNumber.txt','w');
x=fopen('clad_num.txt','w');
fprintf(x,'%d\n',claddingNumber);
fclose(x);

%x=fopen('claddingAbsorption.txt','w');
x=fopen('clad_abs.txt','w');
fprintf(x,'%.50f\n',claddingAbsorption);
fclose(x);

%x=fopen('refractiveIndices.txt','w');
x=fopen('refractive_indices.txt','w');
fprintf(x,'%3.5f\n',refractiveIndices);
fclose(x);

%x=fopen('reflectivities.txt','w');
x=fopen('reflectivities.txt','w');
fprintf(x,'%.50f\n',reflectivities);
fclose(x);

cd(CURRENT_DIR);
end 





%%%%%%%%%%%%%%%%%%%%%%%%% parse_calcPhiASE_output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% takes the output from the CUDA code and fills it into a variable
% assumes that the matrix is saved as a 3D-matrix where the first line 
% denotes the dimensions and the second line is the whole data
%
% @param FOLDER the folder which contains the output files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mseValues,  raysUsedPerSample, phiASE] = parse_calcPhiASE_output (FOLDER)
CURRENT_DIR = pwd;
cd (FOLDER);
fid = fopen('phi_ASE.txt');
arraySize = str2num(fgetl(fid));
phiASE = str2num(fgetl(fid));
phiASE = reshape(phiASE,arraySize);
fclose(fid);

fid=fopen('mse_values.txt');
arraySize=str2num(fgetl(fid));
mseValues = str2num(fgetl(fid));
mseValues = reshape(mseValues,arraySize);
fclose(fid);


fid = fopen('N_rays.txt');
arraySize = str2num(fgetl(fid));
raysUsedPerSample = str2num(fgetl(fid));
raysUsedPerSample = reshape(raysUsedPerSample,arraySize);
fclose(fid);

cd (CURRENT_DIR);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% clean_IO_files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% deletes the temporary folder and the dndt_ASE.txt
% 
% @param TMP_FOLDER the folder to remove
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clean_IO_files (TMP_FOLDER)
A = exist(TMP_FOLDER,'dir'); % continue only if the folder exists!
if A == 7
  s = warning;
  warning off all;

  isOctave = exist('OCTAVE_VERSION') ~= 0;
  if(isOctave)
    confirm_recursive_rmdir (0);
  else
    rmdir(TMP_FOLDER,'s');
  end

  warning(s);
end

end
