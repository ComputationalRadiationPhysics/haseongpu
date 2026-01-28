import vtk
import numpy as np

def vtk_wedge(file_n, phi_ASE, p, t_int, mesh_z, z_mesh):
    # Ensure file extension is .vtk
    if not file_n.endswith('.vtk'):
        file_n += '.vtk'

    num_p_layer = p.shape[0]
    num_t_layer = t_int.shape[0]

    # 1. Create the VTK Points object
    vtk_points = vtk.vtkPoints()
    for i_z in range(mesh_z):
        z = i_z * z_mesh
        for i_p in range(num_p_layer):
            vtk_points.InsertNextPoint(p[i_p, 0], p[i_p, 1], z)

    # 2. Create the Grid and Connectivity
    u_grid = vtk.vtkUnstructuredGrid()
    u_grid.SetPoints(vtk_points)

    for i_z in range(mesh_z - 1):
        layer_offset = i_z * num_p_layer
        next_layer_offset = (i_z + 1) * num_p_layer
        
        for tri in t_int:
            wedge = vtk.vtkWedge()
            # Bottom Triangle
            wedge.GetPointIds().SetId(0, tri[0] + layer_offset)
            wedge.GetPointIds().SetId(1, tri[1] + layer_offset)
            wedge.GetPointIds().SetId(2, tri[2] + layer_offset)
            # Top Triangle
            wedge.GetPointIds().SetId(3, tri[0] + next_layer_offset)
            wedge.GetPointIds().SetId(4, tri[1] + next_layer_offset)
            wedge.GetPointIds().SetId(5, tri[2] + next_layer_offset)
            
            u_grid.InsertNextCell(wedge.GetCellType(), wedge.GetPointIds())

    # 3. Add Scalars (phi_ASE)
    # VTK needs data in a vtkFloatArray
    scalars = vtk.vtkFloatArray()
    scalars.SetName("phi_ASE")
    
    # Flatten the data to match point order (layer by layer)
    flat_data = phi_ASE.flatten(order='F')
    for val in flat_data:
        scalars.InsertNextValue(val)
        
    u_grid.GetPointData().SetScalars(scalars)

    # 4. Write to File
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName(file_n)
    writer.SetInputData(u_grid)
    writer.Write()
    
    print(f"File saved as {file_n}")