import numpy as np
from skimage.graph import route_through_array
image = np.array ([[0, 1, 0, 1, 1, 0],
                   [2, 1, 1, 0, 0, 1],
                   [0, 6, 0, 2, 0, 1],
                   [25, 6, 0, 25, 99, 99],
                   [0, 0, 0, 0, 0, 2],
                   [0, 6, 5, 3, 0, 0]])

indices, weight = route_through_array(image, (0, 0), (5,5))
indices = np.array(indices).T
path = np.zeros_like(image)
path[indices[0], indices[1]] = 1
print path
print weight