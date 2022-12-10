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
