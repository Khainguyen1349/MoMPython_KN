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
