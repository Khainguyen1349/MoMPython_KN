import numpy as np
from scipy.spatial import Delaunay#, ConvexHull


def creatHull(hull_bound_points,mesh_resolution):
    hull_bound = Delaunay(hull_bound_points)
    mesh_points = []
    for i in range(min(i[0] for i in hull_bound_points),max(i[0] for i in hull_bound_points)+1,mesh_resolution):
        for j in range(min(i[1] for i in hull_bound_points),max(i[1] for i in hull_bound_points)+1,mesh_resolution):
            mesh_points.append([i,j])
    mesh_ini = np.array(mesh_points)
    mesh = np.array([p_ for p_ in mesh_ini if hull_bound.find_simplex(p_)>=0]) #check if points are inside hull bound
    tri = Delaunay(mesh)
    return Delaunay_temp(tri.points,tri.simplices)
