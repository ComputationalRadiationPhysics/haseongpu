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
"""
# Load pt.mat file
mat_data = scipy.io.loadmat('pt.mat')
p = mat_data['p']
t = mat_data['t']

set_variables(p, t)  # Call function

print("variable1.mat file generated successfully.")
"""

