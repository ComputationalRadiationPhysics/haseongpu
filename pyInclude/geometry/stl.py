# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

from __future__ import annotations

from pathlib import Path
import struct

import numpy as np


def _read_stl_triangles(path):
    raw = Path(path).read_bytes()
    if len(raw) >= 84:
        tri_count = struct.unpack_from("<I", raw, 80)[0]
        if 84 + tri_count * 50 == len(raw):
            triangles = []
            offset = 84
            for _ in range(tri_count):
                offset += 12
                tri = []
                for _ in range(3):
                    tri.append(struct.unpack_from("<fff", raw, offset))
                    offset += 12
                offset += 2
                triangles.append(tri)
            return triangles

    triangles = []
    current = []
    for line in raw.decode("utf-8", errors="strict").splitlines():
        tokens = line.strip().split()
        if len(tokens) == 4 and tokens[0].lower() == "vertex":
            current.append(tuple(float(v) for v in tokens[1:]))
            if len(current) == 3:
                triangles.append(current)
                current = []
    if not triangles:
        raise ValueError(f"could not parse STL triangles from {path}")
    return triangles


def from_stl(path):
    triangles_3d = _read_stl_triangles(path)
    z = [vertex[2] for tri in triangles_3d for vertex in tri]
    if not np.allclose(z, z[0]):
        raise ValueError("only planar STL meshes are supported for now")

    index = {}
    points = []
    triangles = []
    for tri in triangles_3d:
        ids = []
        for x, y, _ in tri:
            key = (float(x), float(y))
            if key not in index:
                index[key] = len(points)
                points.append(key)
            ids.append(index[key])
        triangles.append(ids)
    return np.asarray(points, dtype=np.float64), np.asarray(triangles, dtype=np.uint32)
