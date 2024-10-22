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


% VTK Wedge (=13) Data writing
% Daniel Albach                                         2009/05/29

function vtk_wedge(file_n, phi_ASE, p, t_int, mesh_z, z_mesh)

% open a file
filename = file_n;

fid = fopen(filename, 'wt');
if fid == -1
    error('Cannot open file for writing.');
end

% newline character
nl = sprintf('\n');

size_p = size(phi_ASE,1);
size_p2 = size(phi_ASE,1)*mesh_z;

fwrite(fid, ['# vtk DataFile Version 2.0' nl 'Wedge example ' nl 'ASCII' nl ...
    'DATASET UNSTRUCTURED_GRID' nl ...
    'POINTS ' num2str(size_p2) ' float' nl]);

% now write the data-coordinates
% create the cells with the z-coodinates at first, then horizontaly
% concatenate
h = zeros(size(p,1),1);
% first level is at zero level
h(:) = z_mesh;

for i_z = 1:mesh_z
    v = horzcat(p,(h.*(i_z-1)));
    for i_v=1:size_p
        fprintf(fid,'%f %f %f', v(i_v,1), v(i_v,2), v(i_v,3));
        fwrite(fid, nl);
    end
end


% now write the cell connections
fwrite(fid, ['CELLS ' num2str((mesh_z-1)*size(t_int,1)) ' ' num2str((mesh_z-1)*7*size(t_int,1) ) nl]);
% in each line: Number of Points points points in c-style indexing
t_0 = int32(zeros(size(t_int,1),1));
t_0(:) = 6;
t_1 = int32(zeros(size(t_int,1),size(t_int,2)));
t_1 = t_int;
t_2 = int32(zeros(size(t_int,1),size(t_int,2)));

for i_z=1:(mesh_z-1)
    
    t_2 = t_1 + size(p,1);
    
    tp = horzcat(t_0, t_1, t_2);
    
    for i_v=1:size(t_int,1)
        fprintf(fid,'%i %i %i %i %i %i %i', tp(i_v,1), tp(i_v,2), tp(i_v,3), tp(i_v,4), tp(i_v,5), tp(i_v,6), tp(i_v,7));
        fwrite(fid, nl);
    end
    
    t_1 = t_2;
end

% now we hve to write the cell types
fwrite(fid, ['CELL_TYPES ' num2str((mesh_z-1)*size(t_int,1)) nl]);
for i_ct=1:(size(t_int,1)*(mesh_z-1))
    fwrite(fid, ['13' nl]);
end

% now write the point data
fwrite(fid, ['POINT_DATA ' num2str(mesh_z*size(p,1)) nl]);
fwrite(fid, ['SCALARS scalars float 1' nl]);
fwrite(fid, ['LOOKUP_TABLE default' nl]);

% now write the PHI_ASE data
for i_z=1:mesh_z
    v = phi_ASE(:,i_z);
    for i_v=1:size_p
        fprintf(fid,'%f', v(i_v,1));
        fwrite(fid, nl);
    end
end

% close the file
fclose(fid);

end