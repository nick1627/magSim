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

# NField.plot3DField(-1.5*Rn, 1.5*Rn, -1.5*Rn, 1.5*Rn, -1.5*Rn, 1.5*Rn, scale = 100000000)
# #NField.plot2DField(np.array([Rn, 0, 0]), np.array([0, Rn, 0]), 6)

#========================================================================================

Ru = 25600000           #radius of Uranus in metres
g = np.array([[11278, 10928, 0], [-9648, -12284, 1453]])
h = np.array([[0, -16049, 0], [0, 6405, 4220]])

UField = SHField(Ru, g, h, 0, 0)

#Plot results of 16 points on the x=0 plane
x, y, z, u, v, w = UField.plot3DField(-1.5*Ru, 1.5*Ru, -1.5*Ru, 1.5*Ru, -1.5*Ru, 1.5*Ru, 4, Ru, planetaryFilter=False)
#x, y, z, u, v, w = UField.plot2DField("x", 0, 16, Ru/4, Ru)

print(x/Ru)
print(y/Ru)
print(z/Ru)
print(u)
print(v)
print(w)

#========================================================================================





plt.show()