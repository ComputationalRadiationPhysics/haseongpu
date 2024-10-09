##########################################################################
 # Copyright 2013 Daniel Albach, Erik Zenker, Carlchristian Eckert
 # 
 #  This file is part of HASEonGPU
 # 
 #  HASEonGPU is free software: you can redistribute it and/or modify
 #  it under the terms of the GNU General Public License as published by
 #  the Free Software Foundation, either version 3 of the License, or
 #  (at your option) any later version.
 # 
 #  HASEonGPU is distributed in the hope that it will be useful,
 #  but WITHOUT ANY WARRANTY; without even the implied warranty of
 #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 #  GNU General Public License for more details.
 # 
 #  You should have received a copy of the GNU General Public License
 #  along with HASEonGPU.
 #  If not, see <http://www.gnu.org/licenses/>.
 ##########################################################################


#  VTK Wedge (=13) Data writing
#  Daniel Albach                                         2009/05/29


import numpy as np

def vtk_wedge(file_n, phi_ASE, p, t_int, mesh_z, z_mesh):
    """
    Write VTK wedge data to a file.

    Parameters:
    file_n : str
        Name of the file to write to.
    phi_ASE : numpy.ndarray
        Array containing ASE data.
    p : numpy.ndarray
        Array containing point coordinates.
    t_int : numpy.ndarray
        Array containing cell connectivity.
    mesh_z : int
        Number of z-levels in the mesh.
    z_mesh : float
        Mesh size in the z-direction.
    """
    
    # Open the file for writing
    with open(file_n, 'wt') as fid:
        nl = '\n'

        size_p = phi_ASE.shape[0] # number of points size p from first dim of phi_ase
        size_p2 = size_p * mesh_z # total number of points in entire mesh
        
        # Write header
        fid.write(f'# vtk DataFile Version 2.0{nl}')
        fid.write('Wedge example{nl}')
        fid.write('ASCII{nl}')
        fid.write('DATASET UNSTRUCTURED_GRID{nl}')
        fid.write(f'POINTS {size_p2} float{nl}')
        
        # Write data-coordinates
        h = np.zeros((p.shape[0], 1)) # h array having same no. of rows and single col
        h[:] = z_mesh # all value in z_mesh represent z coordinate
        
        # loop over each z level, concatenate p point with scaled z corrdinate to create 3D coord v
        # write each point cord in v to the file
        for i_z in range(mesh_z):
            v = np.hstack((p, h * (i_z)))
            for i_v in range(size_p):
                fid.write(f'{v[i_v, 0]} {v[i_v, 1]} {v[i_v, 2]}{nl}')
        
        # Write cell connections
        num_cells = (mesh_z - 1) * t_int.shape[0] # no of z plane -1 * no. of triangle t_int
        num_indices = (mesh_z - 1) * 7 * t_int.shape[0] # total indices, each cell has 7 (6 vertices and 1 for number of vertices)
        fid.write(f'CELLS {num_cells} {num_indices}{nl}')
        
        t_0 = np.zeros((t_int.shape[0], 1), dtype=np.int32) # t_0 has 6 value, points per cell, int32 for cell indexing
        t_0[:] = 6
        t_1 = t_int.astype(np.int32)
        
        # Loop over each z level except last one, t_2 tell cell connectivity incices by number of points p.shape[0]
        # concatenate t_0, t_1, t_2 to get full connectivity tp
        # update t_1 to be t_2 for next iteration
        for i_z in range(mesh_z - 1):
            t_2 = t_1 + p.shape[0]
            tp = np.hstack((t_0, t_1, t_2))
            
            for i_v in range(t_int.shape[0]):
                fid.write(f'{tp[i_v, 0]} {tp[i_v, 1]} {tp[i_v, 2]} {tp[i_v, 3]} {tp[i_v, 4]} {tp[i_v, 5]} {tp[i_v, 6]}{nl}')
            
            t_1 = t_2
        
        # Write cell types 
        ## write the cell type header, which tell number of cell, write cell type to each cell.
        ## 13 corresponds to wedge in VTK
        fid.write(f'CELL_TYPES {num_cells}{nl}')
        for _ in range(num_cells):
            fid.write(f'13{nl}')
        
        # Write point data
        fid.write(f'POINT_DATA {mesh_z * size_p}{nl}')
        fid.write('SCALARS scalars float 1{nl}')
        fid.write('LOOKUP_TABLE default{nl}')
        
        # Write PHI_ASE data
        ## extract the scalar values for the current z-level from phi_ASE
        ## write each scalar value to file
        for i_z in range(mesh_z):
            v = phi_ASE[:, i_z]
            for i_v in range(size_p):
                fid.write(f'{v[i_v]}{nl}')
