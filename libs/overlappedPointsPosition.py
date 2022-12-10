def overlappedPointsPosition(points,points_list):
    olPP = []
    for i in points_list:
        x = (points == i).all(axis=1).nonzero()
        olPP.append(x[0][0])
    return olPP
