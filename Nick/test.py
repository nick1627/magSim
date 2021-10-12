"""
This is the testing script
"""
import numpy as np
import matplotlib.pyplot as plt
from fields import SHField


import fields as f




#Plotting field of Neptune since it's a bit less weird
Rn = 24765000 #radius of neptune in metres
#start with dipole only
g = np.array([[0.09732, 0.03220]])
h = np.array([[0, -0.09889]])

field1 = SHField(Rn, g, h, 0, 0)

#field1.plot3DField(-1.5*Rn, 1.5*Rn, -1.5*Rn, 1.5*Rn, -1.5*Rn, 1.5*Rn, scale = 100000000)
field1.plot2DField(np.array([Rn, 0, 0]), np.array([0, Rn, 0]), 6)

# import matplotlib.pyplot as plt
# import numpy as np

# ax = plt.figure().add_subplot(projection='3d')

# # Make the grid
# x, y, z = np.meshgrid(np.arange(-0.5, 1, 0.5),
#                       np.arange(-0.5, 1, 0.5),
#                       np.arange(-0.5, 1, 0.5))

# print(np.shape(x))

# # Make the direction data for the arrows
# u = np.ones(np.shape(x))
# v = np.ones(np.shape(x))
# w = np.ones(np.shape(x))

# ax.quiver(x, y, z, u, v, w, length=0.1, normalize=True)



plt.show()