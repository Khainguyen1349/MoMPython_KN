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
