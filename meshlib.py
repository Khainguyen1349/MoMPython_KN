## Mesh creation library
import numpy as np
from scipy.spatial import Delaunay#, ConvexHull

class Delaunay_temp:
    def __init__(self,points,simplices):
        self.points = points #(number_of_points,3)
        self.simplices = simplices #(number_of_triangles,3)
    def deletePoint(self,del_point):
        self.points = np.delete(self.points,del_point,0)
        self.simplices = np.delete(self.simplices,np.where(self.simplices == del_point)[0],0)
        for i in range(self.simplices.shape[0]):
            for j in range(3):
                if self.simplices[i,j] >= del_point:
                    self.simplices[i,j] -= 1

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

def overlappedPoints(points, edgecommon):
    oPoints = []
    for i in points:
        if (edgecommon[0][0] == edgecommon[1][0]):
            if ((i[0] == edgecommon[0][0]) and ((i[1] - edgecommon[0][1])*(i[1] - edgecommon[1][1])<=0)):
                oPoints.append([i[0],i[1]])
        elif (edgecommon[0][1] == edgecommon[1][1]):
            if((i[1] == edgecommon[0][1]) and ((i[0] - edgecommon[0][0])*(i[0] - edgecommon[1][0])<=0)):
                oPoints.append([i[0],i[1]])
        elif(((i[0] - edgecommon[0][0])/(edgecommon[1][0]-edgecommon[0][0]) - (i[1] - edgecommon[0][1])/(edgecommon[1][1]-edgecommon[0][1])<10e-3) and ((i[1] - edgecommon[0][1])*(i[1] - edgecommon[1][1])<=0)):
            oPoints.append([i[0],i[1]])
    return np.array(oPoints)

def overlappedPointsPosition(points,points_list):
    olPP = []
    for i in points_list:
        x = (points == i).all(axis=1).nonzero()
        olPP.append(x[0][0])
    return olPP

def extentHull(tri_origine,hull_bound_points_ext,mesh_resolution_ext,commonEdge):
    
    # find overlapped point on tri_original
    olP = overlappedPoints(tri_origine.points, commonEdge)
##    print(olP)
    # find the positions of the ovelapped points
    olPP = overlappedPointsPosition(tri_origine.points,olP)
    # creat hull
    hull_bound_ext = Delaunay(hull_bound_points_ext)
    
    mesh_points = []
    for i in range(min(i[0] for i in hull_bound_points_ext),max(i[0] for i in hull_bound_points_ext)+1,mesh_resolution_ext):
        for j in range(min(i[1] for i in hull_bound_points_ext),max(i[1] for i in hull_bound_points_ext)+1,mesh_resolution_ext):
            mesh_points.append([i,j])
    mesh_ini = np.array(mesh_points)
    
    # find the overlapped point on new tri_
    olP2 = overlappedPoints(mesh_ini, commonEdge)
    # find the positions of the new ovelapped points
    olPP2 = overlappedPointsPosition(mesh_ini,olP2)
    # delete these nex overlapped points
    mesh_ini_del = np.delete(mesh_ini, olPP2, 0)
    # add the ovelapped points in the original tri
    if len(olPP2) > 0:
        mesh_ini = np.append(olP,mesh_ini_del, axis=0)
    # check if inside hull
    mesh_add = np.array([p_ for p_ in mesh_ini if hull_bound_ext.find_simplex(p_)>=0])
    # the extended tri
    tri_ext = Delaunay(mesh_add)
    
    # points in the merged tri
    tri_final_points = np.concatenate((tri_origine.points, mesh_ini_del), axis=0)
    
    tri_ext_simplices = tri_ext.simplices + tri_origine.points.shape[0] - len(olPP)
    for i in range(len(olPP)):
        pos = np.where(tri_ext_simplices == (i + tri_origine.points.shape[0] - len(olPP)))
        for j in range(len(pos[0])):
            tri_ext_simplices[pos[0][j]][pos[1][j]] = olPP[i]
            
    tri_final_simplices = np.concatenate((tri_origine.simplices, tri_ext_simplices), axis=0)
    
    return Delaunay_temp(tri_final_points,tri_final_simplices)

def closeLoopHull(tri_origine,hull_bound_points_ext,mesh_resolution_ext,commonEdge1,commonEdge2):
    
    # find overlapped point on tri_original
    olP = np.concatenate((overlappedPoints(tri_origine.points, commonEdge1), overlappedPoints(tri_origine.points, commonEdge2)), axis=0)
    # find the positions of the ovelapped points
    olPP = overlappedPointsPosition(tri_origine.points,olP)
    # creat hull
    hull_bound_ext = Delaunay(hull_bound_points_ext)
    
    mesh_points = []
    for i in range(min(i[0] for i in hull_bound_points_ext),max(i[0] for i in hull_bound_points_ext)+1,mesh_resolution_ext):
        for j in range(min(i[1] for i in hull_bound_points_ext),max(i[1] for i in hull_bound_points_ext)+1,mesh_resolution_ext):
            mesh_points.append([i,j])
    mesh_ini = np.array(mesh_points)

    # find the overlapped point on new tri_
    olP2 = np.concatenate((overlappedPoints(mesh_ini, commonEdge1), overlappedPoints(mesh_ini, commonEdge2)), axis=0)
    # find the positions of the new ovelapped points
    olPP2 = overlappedPointsPosition(mesh_ini,olP2)
    # delete these nex overlapped points
    mesh_ini_del = np.delete(mesh_ini, olPP2, 0)
    # add the ovelapped points in the original tri
    if len(olPP2) > 0:
        mesh_ini = np.append(olP,mesh_ini_del, axis=0)
    # check if inside hull
    mesh_add = np.array([p_ for p_ in mesh_ini if hull_bound_ext.find_simplex(p_)>=0])
    # the extended tri
    tri_ext = Delaunay(mesh_add)
    
    # points in the merged tri
    tri_final_points = np.concatenate((tri_origine.points, mesh_ini_del), axis=0)
    
    tri_ext_simplices = tri_ext.simplices + tri_origine.points.shape[0] - len(olPP)
    for i in range(len(olPP)):
        pos = np.where(tri_ext_simplices == (i + tri_origine.points.shape[0] - len(olPP)))
        for j in range(len(pos[0])):
            tri_ext_simplices[pos[0][j]][pos[1][j]] = olPP[i]
            
    tri_final_simplices = np.concatenate((tri_origine.simplices, tri_ext_simplices), axis=0)
    
    return Delaunay_temp(tri_final_points,tri_final_simplices)

def addSeparatedShapes(tri1, tri2):
    tri_final_points = np.concatenate((tri1.points, tri2.points), axis=0)
    #print("size of simplices",tri1.simplices.shape)
    tri_final_simplices = np.concatenate((tri1.simplices, tri2.simplices + tri1.points.shape[0]), axis=0)
    return Delaunay_temp(tri_final_points,tri_final_simplices)

def mirror2DShapes(tri,axis,line):
    ## only accept that the original shape on the positive side of the mirroring axis
    if axis == "X":
        if min(tri.points[0,:]) > line:
            return addSeparatedShapes(tri, Delaunay_temp(tri.points*np.array([[-1],[1]]) + np.array([[2*line],[0]]),tri.simplices))
        elif min(tri.points[0,:]) == line:
            print("Error!")
            return tri
    else:
        if min(tri.points[1,:]) > line:
            return addSeparatedShapes(tri, Delaunay_temp(tri.points*np.array([[1],[-1]]) + np.array([[0],[2*line]]),tri.simplices))
        elif min(tri.points[1,:]) == line:
            print("Error!")
            return tri

def to3D(tri2D,scaleX,scaleY,rotX,rotY,rotZ,dX,dY,dZ):
    #rotX,rotY,rotZ in radiant
    translate = np.array([[dX,dY,dZ]])
    scale = np.array([[scaleX,scaleY]]).T
    
    RotX = np.array([[1,0,0],
                     [0,np.cos(rotX),-np.sin(rotX)],
                     [0,np.sin(rotX),np.cos(rotX)]])
    RotY = np.array([[np.cos(rotY),0,np.sin(rotY)],
                     [0,1,0],
                     [-np.sin(rotY),0,np.cos(rotY)]])
    RotZ = np.array([[np.cos(rotZ),-np.sin(rotZ),0],
                     [np.sin(rotZ),np.cos(rotZ),0],
                     [0,0,1]])
    perm = np.matmul(np.matmul(RotX,RotY),RotZ)

    tri3D = tri2D
    points = np.zeros((3,tri3D.points.shape[0]))
    points[0:2,:] = tri3D.points.T*scale
    tri3D.points = (np.matmul(points.T,perm) + translate)
    return tri3D


    
