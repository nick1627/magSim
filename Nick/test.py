"""
This is the testing script
"""
import sys, os
sys.path.insert(0, os.getcwd())

import numpy as np
import matplotlib.pyplot as plt
from fields import SHField
import tools



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


tools.saveBField(x, y, z, u, v, w, "Output/quadrupole_nick3.npz")

a, b, c, d, e, f = tools.loadBField("Output/quadrupole_nick3.npz")

print(x - a)

# print(x/Ru)
# print(y/Ru)
# print(z/Ru)
# print(u)
# print(v)
# print(w)











#========================================================================================

# Ru = 25600000           #radius of Uranus in metres
# g = np.array([[0, 0, 0], [1, 0, 0]])
# h = np.array([[0, 0, 0], [0, 0, 0]])



# UField = SHField(Ru, g, h, 0, 0)

# #Plot results of 16 points on the x=0 plane
# #x, y, z, u, v, w = UField.plot3DField(-1.5*Ru, 1.5*Ru, -1.5*Ru, 1.5*Ru, -1.5*Ru, 1.5*Ru, 4, Ru, planetaryFilter=False)
# x, y, z, u, v, w = UField.plot2DField("x", 0, 20, Ru*(3/19), Ru)

# #print(x/Ru)
# # print(y/Ru)
# # print(z/Ru)
# # print(u)
# # print(v)
# # print(w)



#plt.show()