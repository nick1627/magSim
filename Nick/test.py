"""
This is the testing script
"""
import numpy as np
import matplotlib.pyplot as plt
from fields import SHField


import fields as f



#========================================================================================
# #Plotting field of Neptune since it's a bit less weird
# Rn = 24765000 #radius of neptune in metres
# #start with dipole only
# g = np.array([[0.09732, 0.03220]])
# h = np.array([[0, -0.09889]])

# NField = SHField(Rn, g, h, 0, 0)

# #field1.plot3DField(-1.5*Rn, 1.5*Rn, -1.5*Rn, 1.5*Rn, -1.5*Rn, 1.5*Rn, scale = 100000000)
# NField.plot2DField(np.array([Rn, 0, 0]), np.array([0, Rn, 0]), 6)

#========================================================================================

Ru = 25600000           #radius of Uranus in metres
g = np.array([[11893, 11579]])
h = np.array([[0, -15684]])

UField = SHField(Ru, g, h, 0, 0)

#Plot results of 16 points on the x=0 plane
x, y, z, u, v, w = UField.plot2DField(np.array([0, Ru, 0]), np.array([0, 0, Ru]), 4, scale = 1) 

print(y/Ru)
print(z/Ru)
print(u)
print(v)
print(w)

#========================================================================================





plt.show()