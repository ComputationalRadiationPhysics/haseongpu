####################################################################################
# Copyright 2013 Daniel Albach, Erik Zenker, Carlchristian Eckert
#
# This file is part of HASEonGPU
#
# HASEonGPU is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HASEonGPU is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HASEonGPU.
# If not, see <http://www.gnu.org/licenses/>.
###################################################################################

import numpy as np
import scipy
from scipy.io import savemat
"""
def set_variables(p, t):
    t = t - 1  # Convert to zero-based indexing for Python
    length_t = len(t)

    # Step 1: Calculate the index of neighboring cells
    sorted_int = np.zeros((length_t, 3), dtype=int)

    # Now search the cells which have two of the three points in common
    for i_point in range(length_t):
        point_1, point_2, point_3 = t[i_point]
        for i_test in range(length_t):
            if i_test != i_point:
                test_points = t[i_test]
                if len(set([point_1, point_2]) & set(test_points)) == 2:
                    sorted_int[i_point, 0] = i_test + 1
                if len(set([point_1, point_3]) & set(test_points)) == 2:
                    sorted_int[i_point, 1] = i_test + 1
                if len(set([point_2, point_3]) & set(test_points)) == 2:
                    sorted_int[i_point, 2] = i_test + 1

    # Step 2: Make a list of the forbidden faces of the triangles
    forbidden = np.zeros((length_t, 3), dtype=int)
    for i_point in range(length_t):
        face_1 = sorted_int[i_point, 0] - 1
        if face_1 != -1:
            if sorted_int[face_1, 0] == i_point + 1:
                forbidden[i_point, 0] = 1
            elif sorted_int[face_1, 1] == i_point + 1:
                forbidden[i_point, 0] = 2
            elif sorted_int[face_1, 2] == i_point + 1:
                forbidden[i_point, 0] = 3

        face_2 = sorted_int[i_point, 1] - 1
        if face_2 != -1:
            if sorted_int[face_2, 0] == i_point + 1:
                forbidden[i_point, 1] = 1
            elif sorted_int[face_2, 1] == i_point + 1:
                forbidden[i_point, 1] = 2
            elif sorted_int[face_2, 2] == i_point + 1:
                forbidden[i_point, 1] = 3

        face_3 = sorted_int[i_point, 2] - 1
        if face_3 != -1:
            if sorted_int[face_3, 0] == i_point + 1:
                forbidden[i_point, 2] = 1
            elif sorted_int[face_3, 1] == i_point + 1:
                forbidden[i_point, 2] = 2
            elif sorted_int[face_3, 2] == i_point + 1:
                forbidden[i_point, 2] = 3

    # Step 3: Calculate the centers of the triangles via barycentric coordinates
    center_vec = np.array([1 / 3, 1 / 3, 1 / 3], dtype=np.float32)
    #x_center = np.zeros(length_t)
    #y_center = np.zeros(length_t)
    x_center = np.zeros((length_t, 1), dtype=np.float64)  # Column vector (vertical)
    y_center = np.zeros((length_t, 1), dtype=np.float64)  # Column vector (vertical)

    for i_points in range(length_t):
        Px = [p[t[i_points, 0], 0], p[t[i_points, 1], 0], p[t[i_points, 2], 0]]
        Py = [p[t[i_points, 0], 1], p[t[i_points, 1], 1], p[t[i_points, 2], 1]]

        # Calculate the center using barycentric coordinates
        #x_center[i_points] = np.dot(Px, center_vec)
        #y_center[i_points] = np.dot(Py, center_vec)
        x_center[i_points, 0] = np.dot(Px, center_vec)
        y_center[i_points, 0] = np.dot(Py, center_vec)

    # Step 4: Calculate the normal vectors of the different combinations
    normals_x = np.zeros((length_t, 3))
    normals_y = np.zeros((length_t, 3))
    normals_z = np.zeros((length_t, 3))
    normals_p = np.zeros((length_t, 3), dtype=int)

    for i_points in range(length_t):
        # Normals to points 1-2
        vec_1 = np.array([p[t[i_points, 0], 0], p[t[i_points, 0], 1], 0]) - np.array([p[t[i_points, 1], 0], p[t[i_points, 1], 1], 0])
        vec_2 = np.array([p[t[i_points, 0], 0], p[t[i_points, 0], 1], 0]) - np.array([p[t[i_points, 0], 0], p[t[i_points, 0], 1], 0.1])
        vec_cross = np.cross(vec_1, vec_2)
        vec_cross /= np.linalg.norm(vec_cross)
        normals_x[i_points, 0] = vec_cross[0]
        normals_y[i_points, 0] = vec_cross[1]
        normals_z[i_points, 0] = vec_cross[2]
        normals_p[i_points, 0] = t[i_points, 0]

        # Normals to points 1-3
        vec_1 = np.array([p[t[i_points, 0], 0], p[t[i_points, 0], 1], 0]) - np.array([p[t[i_points, 2], 0], p[t[i_points, 2], 1], 0])
        vec_cross = np.cross(vec_1, vec_2)
        vec_cross /= np.linalg.norm(vec_cross)
        normals_x[i_points, 1] = vec_cross[0]
        normals_y[i_points, 1] = vec_cross[1]
        normals_z[i_points, 1] = vec_cross[2]
        normals_p[i_points, 1] = t[i_points, 0]

        # Normals to points 2-3
        vec_1 = np.array([p[t[i_points, 1], 0] - p[t[i_points, 2], 0], p[t[i_points, 1], 1] - p[t[i_points, 2], 1], 0])
        vec_cross = np.cross(vec_1, vec_2)
        vec_cross /= np.linalg.norm(vec_cross)
        normals_x[i_points, 2] = vec_cross[0]
        normals_y[i_points, 2] = vec_cross[1]
        normals_z[i_points, 2] = vec_cross[2]
        normals_p[i_points, 2] = t[i_points, 1]

    # Step 5: Calculate the surface of the triangles
    surface = np.zeros((len(t), 1))
    for i_point in range(len(t)):
        a = ((p[t[i_point, 0], 0] - p[t[i_point, 1], 0]) ** 2 + (p[t[i_point, 0], 1] - p[t[i_point, 1], 1]) ** 2)
        b = ((p[t[i_point, 2], 0] - p[t[i_point, 1], 0]) ** 2 + (p[t[i_point, 2], 1] - p[t[i_point, 1], 1]) ** 2)
        c = ((p[t[i_point, 0], 0] - p[t[i_point, 2], 0]) ** 2 + (p[t[i_point, 0], 1] - p[t[i_point, 2], 1]) ** 2)
        surface[i_point] = np.sqrt(1 / 16 * (4 * a * c - (a + c - b) ** 2))

    # Convert variables to int32 as necessary for compatibility
    t_int = t.astype(np.int32)
    sorted_int = sorted_int.astype(np.int32) - 1
    normals_p = normals_p.astype(np.int32)
    forbidden = forbidden.astype(np.int32) -1

    # Save the variables to a MAT file
    savefile = 'variable.mat'
    savemat(savefile, {
        'normals_x': normals_x,
        'normals_y': normals_y,
        'normals_z': normals_z,
        'sorted_int': sorted_int,
        't_int': t_int,
        'x_center': x_center,
        'y_center': y_center,
        'surface': surface,
        'normals_p': normals_p,
        'forbidden': forbidden
    })

# Load pt.mat file
mat_data = scipy.io.loadmat('pt.mat')
p = mat_data['p']
t = mat_data['t']

set_variables(p, t)  # Call function

print("variable1.mat file generated successfully.")
"""

####################################################################################
# Copyright 2013 Daniel Albach, Erik Zenker, Carlchristian Eckert
#
# This file is part of HASEonGPU
#
# HASEonGPU is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HASEonGPU is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HASEonGPU.
# If not, see <http://www.gnu.org/licenses/>.
###################################################################################
"""
import numpy as np
from scipy.io import savemat

def set_variables(p, t):
    
    length_t = t.shape[0]

    # Step 1: Calculate neighbor indices
    sorted_int = np.zeros((length_t, 3), dtype=int)
    for i in range(length_t):
        point_1, point_2, point_3 = t[i,:]
        for j in range(length_t):
            if j == i:
                continue
            test_points = t[j,:]
            if len(set([point_1, point_2]) & set(test_points)) == 2:
                sorted_int[i,0] = j + 1  # 1-based temporarily
            if len(set([point_1, point_3]) & set(test_points)) == 2:
                sorted_int[i,1] = j + 1
            if len(set([point_2, point_3]) & set(test_points)) == 2:
                sorted_int[i,2] = j + 1

    # Step 2: Forbidden faces
    forbidden = np.zeros((length_t, 3), dtype=int)
    for i in range(length_t):
        for k in range(3):
            face = sorted_int[i,k] - 1  # convert temporarily to 0-based
            if face < 0:
                forbidden[i,k] = -1
                continue
            for neighbor_face in range(3):
                if sorted_int[face,neighbor_face] - 1 == i:
                    forbidden[i,k] = neighbor_face
                    break

    # Step 3: Triangle centers (barycentric)
    center_vec = np.array([1/3,1/3,1/3], dtype=np.float64)
    x_center = np.zeros(length_t, dtype=np.float64)
    y_center = np.zeros(length_t, dtype=np.float64)
    for i in range(length_t):
        Px = p[t[i,:]-1,0]  # convert to 0-based temporarily
        Py = p[t[i,:]-1,1]
        x_center[i] = np.dot(Px, center_vec)
        y_center[i] = np.dot(Py, center_vec)

    # Step 4: Normals
    normals_x = np.zeros((length_t,3))
    normals_y = np.zeros((length_t,3))
    normals_z = np.zeros((length_t,3))
    normals_p = np.zeros((length_t,3), dtype=int)

    for i in range(length_t):
        pts = t[i,:]-1  # 0-based
        # 1-2
        vec_1 = np.array([p[pts[0],0]-p[pts[1],0], p[pts[0],1]-p[pts[1],1], 0])
        vec_2 = np.array([0,0,0.1])
        cross = np.cross(vec_1, vec_2)
        cross /= np.linalg.norm(cross)
        normals_x[i,0] = cross[0]; normals_y[i,0] = cross[1]; normals_z[i,0] = cross[2]
        normals_p[i,0] = pts[0]

        # 1-3
        vec_1 = np.array([p[pts[0],0]-p[pts[2],0], p[pts[0],1]-p[pts[2],1], 0])
        cross = np.cross(vec_1, vec_2)
        cross /= np.linalg.norm(cross)
        normals_x[i,1] = cross[0]; normals_y[i,1] = cross[1]; normals_z[i,1] = cross[2]
        normals_p[i,1] = pts[0]

        # 2-3
        vec_1 = np.array([p[pts[1],0]-p[pts[2],0], p[pts[1],1]-p[pts[2],1], 0])
        cross = np.cross(vec_1, vec_2)
        cross /= np.linalg.norm(cross)
        normals_x[i,2] = cross[0]; normals_y[i,2] = cross[1]; normals_z[i,2] = cross[2]
        normals_p[i,2] = pts[1]

    # Step 5: Triangle surfaces
    surface = np.zeros(length_t)
    for i in range(length_t):
        pts = t[i,:]-1
        a = (p[pts[0],0]-p[pts[1],0])**2 + (p[pts[0],1]-p[pts[1],1])**2
        b = (p[pts[2],0]-p[pts[1],0])**2 + (p[pts[2],1]-p[pts[1],1])**2
        c = (p[pts[0],0]-p[pts[2],0])**2 + (p[pts[0],1]-p[pts[2],1])**2
        surface[i] = np.sqrt(1/16*(4*a*c-(a+c-b)**2))

    # Step 6: Convert to int32 and zero-based indexing for CUDA
    t_int = t.astype(np.int32) - 1
    sorted_int = sorted_int.astype(np.int32) - 1
    normals_p = normals_p.astype(np.int32)
    forbidden = forbidden.astype(np.int32)  # already zero-based with -1 for no neighbor

    # Save MAT file
    savefile = 'variable.mat'
    savemat(savefile, {
        'normals_x': normals_x,
        'normals_y': normals_y,
        'normals_z': normals_z,
        'sorted_int': sorted_int,
        't_int': t_int,
        'x_center': x_center,
        'y_center': y_center,
        'surface': surface,
        'normals_p': normals_p,
        'forbidden': forbidden
    })

    print("variable.mat file generated successfully.")
"""

####################################################################################
# set_variable.py (HASEonGPU)
# - Compacts/remaps points so all points are referenced by triangles
# - Writes parser-safe arrays:
#     triangleNormalPoint: unsigned int in [0 .. numberOfPoints-1] with max == numberOfPoints-1
#     triangleNeighbors:   int in [-1 .. numberOfTriangles-1] with exact min/max
#     forbiddenEdge:       int in [-1 .. 2] with exact min/max
####################################################################################

import numpy as np
import scipy.io
from scipy.io import savemat


# def compact_mesh(p, t):
#     """
#     p: (N,2) float
#     t: (M,3) int, 1-based (MATLAB style)
#     Returns:
#       p2: compacted points (K,2)
#       t2: remapped triangles (M,3) still 1-based
#     """
#     p = np.asarray(p, dtype=np.float64)
#     t = np.asarray(t, dtype=np.int64)
#
#     # Flatten triangle indices, unique used vertices (still 1-based)
#     used = np.unique(t.reshape(-1))
#     used_sorted = np.sort(used)
#
#     # Build mapping old_index(1-based) -> new_index(1-based)
#     # Use dict for clarity; could be vectorized but this is robust.
#     mapping = {old: new for new, old in enumerate(used_sorted, start=1)}
#
#     # Remap triangles
#     t2 = np.vectorize(mapping.get)(t).astype(np.int32)
#
#     # Subset points (convert used indices to 0-based for slicing)
#     p2 = p[used_sorted - 1, :].astype(np.float64)
#
#     return p2, t2


def set_variables(p, t):
    """
    Builds variable.mat in a way that matches parser.cu assertions.
    """
    # 1) Compact mesh so numberOfPoints matches max index used by triangles
    ## p, t = compact_mesh(p, t) do not compact the mesh this will cause index space problems later on

    # Sanitize
    p = np.asarray(p, dtype=np.float64)
    t = np.asarray(t, dtype=np.int32)

    num_tri = t.shape[0]
    num_pts = p.shape[0]

    # DO not substract -1 as t is already 0 based
    # t0 = t - 1  # shape (num_tri, 3), values 0..num_pts-1

    # =========================================================
    # STEP 1: triangleNeighbors (sorted_int)
    # MATLAB: neighbors stored 1..num_tri, 0 means "none"
    # We will finally store as int32 with -1 meaning none, else 0..num_tri-1
    # =========================================================
    sorted_int = np.zeros((num_tri, 3), dtype=np.int32)  # 0 means "no neighbor" (MATLAB style)

    for i in range(num_tri):
        a, b, c = t[i]
        for j in range(num_tri):
            if i == j:
                continue
            tj = t[j]
            # Share edge (a,b)
            if len({a, b}.intersection(tj)) == 2:
                sorted_int[i, 0] = j + 1
            # Share edge (a,c)
            if len({a, c}.intersection(tj)) == 2:
                sorted_int[i, 1] = j + 1
            # Share edge (b,c)
            if len({b, c}.intersection(tj)) == 2:
                sorted_int[i, 2] = j + 1

    # Convert neighbors to int32 with -1 for boundary, else 0..num_tri-1
    triangleNeighbors = (sorted_int - 1).astype(np.int32)  # 0->-1, (j+1)->j

    # =========================================================
    # STEP 2: forbiddenEdge
    # MATLAB: forbidden in {0,1,2,3} where 0 = none
    # parser.cu wants int range exactly [-1..2] (equals==true)
    # We'll store -1 for none, else 0..2
    # =========================================================
    forbidden = np.zeros((num_tri, 3), dtype=np.int32)  # MATLAB style: 0 means none, else 1..3

    for i in range(num_tri):
        for k in range(3):
            face = sorted_int[i, k] - 1  # neighbor triangle index (0-based) or -1
            if face < 0:
                continue
            # find which edge index (1..3) in neighbor points back to i
            for kk in range(3):
                if (sorted_int[face, kk] - 1) == i:
                    forbidden[i, k] = kk + 1
                    break

    forbiddenEdge = (forbidden - 1).astype(np.int32)  # 0->-1, 1..3 -> 0..2

    # =========================================================
    # STEP 3: triangle centers
    # =========================================================
    x_center = np.zeros(num_tri, dtype=np.float64)
    y_center = np.zeros(num_tri, dtype=np.float64)
    for i in range(num_tri):
        pts = t[i]
        x_center[i] = np.mean(p[pts, 0])
        y_center[i] = np.mean(p[pts, 1])

    # =========================================================
    # STEP 4: normals + triangleNormalPoint
    # parser.cu expects triangleNormalPoint unsigned with min=0 max=num_pts-1 (equals true)
    #
    # We'll follow the original logic:
    #  edge 1-2 uses point 1
    #  edge 1-3 uses point 1
    #  edge 2-3 uses point 2
    #
    # IMPORTANT: after compaction, max point index WILL be used somewhere (because points are only used ones)
    # =========================================================
    normals_x = np.zeros((num_tri, 3), dtype=np.float64)
    normals_y = np.zeros((num_tri, 3), dtype=np.float64)
    normals_z = np.zeros((num_tri, 3), dtype=np.float64)

    triangleNormalPoint = np.zeros((num_tri, 3), dtype=np.uint32)

    vec2 = np.array([0.0, 0.0, 0.1], dtype=np.float64)

    for i in range(num_tri):
        pts = t[i]

        # Edge 1-2
        vec1 = np.array([p[pts[0], 0], p[pts[0], 1], 0.0], dtype=np.float64) - \
               np.array([p[pts[1], 0], p[pts[1], 1], 0.0], dtype=np.float64)
        n = np.cross(vec1, vec2)
        n /= np.linalg.norm(n)
        normals_x[i, 0], normals_y[i, 0], normals_z[i, 0] = n
        triangleNormalPoint[i, 0] = np.uint32(pts[0])

        # Edge 1-3
        vec1 = np.array([p[pts[0], 0], p[pts[0], 1], 0.0], dtype=np.float64) - \
               np.array([p[pts[2], 0], p[pts[2], 1], 0.0], dtype=np.float64)
        n = np.cross(vec1, vec2)
        n /= np.linalg.norm(n)
        normals_x[i, 1], normals_y[i, 1], normals_z[i, 1] = n
        triangleNormalPoint[i, 1] = np.uint32(pts[0])

        # Edge 2-3
        vec1 = np.array([p[pts[1], 0], p[pts[1], 1], 0.0], dtype=np.float64) - \
               np.array([p[pts[2], 0], p[pts[2], 1], 0.0], dtype=np.float64)
        n = np.cross(vec1, vec2)
        n /= np.linalg.norm(n)
        normals_x[i, 2], normals_y[i, 2], normals_z[i, 2] = n
        triangleNormalPoint[i, 2] = np.uint32(pts[1])

    # =========================================================
    # STEP 5: triangle surface (MATLAB formula)
    # =========================================================
    surface = np.zeros(num_tri, dtype=np.float64)
    for i in range(num_tri):
        a = np.sum((p[t[i, 0]] - p[t[i, 1]]) ** 2)
        b = np.sum((p[t[i, 2]] - p[t[i, 1]]) ** 2)
        c = np.sum((p[t[i, 0]] - p[t[i, 2]]) ** 2)
        surface[i] = np.sqrt(1.0 / 16.0 * (4.0 * a * c - (a + c - b) ** 2))

    # =========================================================
    # Save MAT
    # Notes:
    # - Keep triangleNeighbors / forbiddenEdge SIGNED int32 because parser expects -1 allowed.
    # - triangleNormalPoint must be uint32 0..num_pts-1 with exact max == num_pts-1.
    # - trianglePointIndices ("t_int") should be uint32 0..num_pts-1.
    # =========================================================
    savemat("variable.mat", {
        "t_int": t.astype(np.uint32),
        "sorted_int": triangleNeighbors.astype(np.int32),
        "forbidden": forbiddenEdge.astype(np.int32),
        "normals_p": triangleNormalPoint.astype(np.uint32),

        "normals_x": normals_x,
        "normals_y": normals_y,
        "normals_z": normals_z,
        "x_center": x_center,
        "y_center": y_center,
        "surface": surface,
    })

    # Quick diagnostics (helps confirm parser assertions)
    tp = triangleNormalPoint.reshape(-1)
    tn = triangleNeighbors.reshape(-1)
    fe = forbiddenEdge.reshape(-1)

    print("variable.mat generated (parser-safe).")
    print(f"[check] numberOfPoints={num_pts}  numberOfTriangles={num_tri}")
    print(f"[check] triangleNormalPoint min={tp.min()} max={tp.max()} (expect 0..{num_pts-1})")
    print(f"[check] triangleNeighbors   min={tn.min()} max={tn.max()} (expect -1..{num_tri-1})")
    print(f"[check] forbiddenEdge       min={fe.min()} max={fe.max()} (expect -1..2)")


if __name__ == "__main__":
    mat = scipy.io.loadmat("pt.mat")
    set_variables(mat["p"], mat["t"])
