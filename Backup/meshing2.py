import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

points = []
for i in range(9):
    for j in range(9):
        if (i-4)**2 + (j-4)**2 <= 4**2:
            points.append([i - 4,j - 4])


parray = np.array(points)

tri = Delaunay(parray)
plt.triplot(parray[:,0], parray[:,1], tri.simplices)
plt.show()
