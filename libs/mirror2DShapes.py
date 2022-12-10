def mirror2DShapes(tri,axis,line,mirrorOverO):
    ## Strictly required that the original shape on the positive side of the mirroring axis and have at least 2 range of mesh at the mirrored edge
    if axis == "X":
        print("Min X:",min(tri.points[:,0]))
        if min(tri.points[:,0]) > line:
            if mirrorOverO:
                return addSeparatedShapes(tri, Delaunay_temp(tri.points*np.array([-1,-1]) + np.array([2*line,0]),tri.simplices))
            else:
                return addSeparatedShapes(tri, Delaunay_temp(tri.points*np.array([-1,1]) + np.array([2*line,0]),tri.simplices))
        elif min(tri.points[:,0]) == line:
            mirrorEdge = np.array([[line,int(np.min(tri.points[np.where(tri.points[:,0] == line),1]))],
                                   [line,int(np.max(tri.points[np.where(tri.points[:,0] == line),1]))]])
            overlappedPoints(tri.points, mirrorEdge)
            pp_list = overlappedPointsPosition(tri.points,overlappedPoints(tri.points, mirrorEdge))
            if mirrorOverO:
                tri_temp = Delaunay_temp(tri.points*np.array([-1,-1]) + np.array([2*line,0]), tri.simplices)
            else:
                tri_temp = Delaunay_temp(tri.points*np.array([-1,1]) + np.array([2*line,0]), tri.simplices)
            for i in range(len(pp_list)):
                tri_temp.deletePoint(pp_list[i]-i)
                
            tri = addSeparatedShapes(tri, tri_temp)

            hull_bound_points_ext = np.array([[int(max(tri_temp.points[:,0])),mirrorEdge[0,1]],
                                              [mirrorEdge[0,0],mirrorEdge[0,1]],
                                              [int(max(tri_temp.points[:,0])),mirrorEdge[1,1]],
                                               [mirrorEdge[1,0],mirrorEdge[1,1]]])
            #print(hull_bound_points_ext)
            mesh_resolution_ext = int(min([np.abs(tri.points[pp_list[0],1] - tri.points[i,1]) for i in pp_list[1:]]))
            #print(mesh_resolution_ext)
            commonEdge1 = mirrorEdge
            #print(commonEdge1)
            commonEdge2 = np.array([[int(max(tri_temp.points[:,0])),mirrorEdge[0,1]],
                                    [int(max(tri_temp.points[:,0])),mirrorEdge[1,1]]])
            #print(commonEdge2)
            tri = closeLoopHull(tri,hull_bound_points_ext,mesh_resolution_ext,commonEdge1,commonEdge2)
            return tri
        else:
            print("Error: overlapping mirroring is not support!")
    else:
        print("Min Y:",min(tri.points[:,1]))
        if min(tri.points[:,1]) > line:
            if mirrorOverO:
                return addSeparatedShapes(tri, Delaunay_temp(tri.points*np.array([-1,-1]) + np.array([0,2*line]),tri.simplices))
            else:
                return addSeparatedShapes(tri, Delaunay_temp(tri.points*np.array([1,-1]) + np.array([0,2*line]),tri.simplices))
        elif min(tri.points[:,1]) == line:
            mirrorEdge = np.array([[int(np.min(tri.points[np.where(tri.points[:,1] == line),0])),line],
                                   [int(np.max(tri.points[np.where(tri.points[:,1] == line),0])),line]])
            overlappedPoints(tri.points, mirrorEdge)
            pp_list = overlappedPointsPosition(tri.points,overlappedPoints(tri.points, mirrorEdge))
            if mirrorOverO:
                tri_temp = Delaunay_temp(tri.points*np.array([-1,-1]) + np.array([0,2*line]), tri.simplices)
            else:
                tri_temp = Delaunay_temp(tri.points*np.array([1,-1]) + np.array([0,2*line]), tri.simplices)
            for i in range(len(pp_list)):
                tri_temp.deletePoint(pp_list[i]-i)
                
            tri = addSeparatedShapes(tri, tri_temp)

            hull_bound_points_ext = np.array([[mirrorEdge[0,0],int(max(tri_temp.points[:,1]))],
                                              [mirrorEdge[1,0],int(max(tri_temp.points[:,1]))],
                                              [mirrorEdge[0,0],mirrorEdge[0,1]],
                                               [mirrorEdge[1,0],mirrorEdge[1,1]]])
            #print(hull_bound_points_ext)
            mesh_resolution_ext = int(min([np.abs(tri.points[pp_list[0],0] - tri.points[i,0]) for i in pp_list[1:]]))
            #print(mesh_resolution_ext)
            commonEdge1 = mirrorEdge
            #print(commonEdge1)
            commonEdge2 = np.array([[mirrorEdge[0,0],int(max(tri_temp.points[:,1]))],
                                    [mirrorEdge[1,0],int(max(tri_temp.points[:,1]))]])
            #print(commonEdge2)
            tri = closeLoopHull(tri,hull_bound_points_ext,mesh_resolution_ext,commonEdge1,commonEdge2)
            return tri
        else:
            print("Error: overlapping mirroring is not support!")
