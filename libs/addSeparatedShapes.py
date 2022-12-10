def addSeparatedShapes(tri1, tri2):
    tri_final_points = np.concatenate((tri1.points, tri2.points), axis=0)
    #print("size of simplices",tri1.simplices.shape)
    tri_final_simplices = np.concatenate((tri1.simplices, tri2.simplices + tri1.points.shape[0]), axis=0)
    return Delaunay_temp(tri_final_points,tri_final_simplices)
