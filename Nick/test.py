"""
This is the testing script
"""
import sys, os
sys.path.insert(0, os.getcwd())

import numpy as np
import matplotlib.pyplot as plt
from fields import SHField, OTDField
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

# g = np.array([[11278, 10928]])
# h = np.array([[0, -16049]])

# g = np.array([[0, 0]])
# h = np.array([[0, 1]])

g = np.array([[11278, 10928, 0], [-9648, -12284, 1453]])
h = np.array([[0, -16049, 0], [0, 6405, 4220]])

UField = SHField(Ru, g, h, 0, 0)
# UField_dipole = OTDField(np.array([0, 1, 0]))
UField.rotate("Field")

# print(UField.getField(np.array([0, 0, 0])))
#Plot results of 16 points on the x=0 plane
# x, y, z, u, v, w = UField.plot3DField(-1.5*Ru, 1.5*Ru, -1.5*Ru, 1.5*Ru, -1.5*Ru, 1.5*Ru, 4, Ru, planetaryFilter=False)
# a = UField.plot2DField("y", 0, 16, Ru/4, Ru)
# b = UField.plot2DField("y", 0, 16, Ru/4, Ru)
# c = UField.plot2DField("z", 0, 16, Ru/4, Ru)

# print(u[0, 0, 0] - 1.79626500e+03)
# tools.saveBField(x, y, z, u, v, w, "Output/dipole_nick_fieldaligned.npz")

# print(np.sqrt(a[3][4, 4]**2 + a[4][4, 4]**2 + a[5][4, 4]**2))
# print(np.sqrt(b[3][4, 4]**2 + b[4][4, 4]**2 + b[5][4, 4]**2))

# UField.plotDeviationData(5, 5, 10)
# UField.plotLongitudePlanesB(5, 5, 10)

# true_nMax = UField.nMax
# UField.nMax = 1

# Bd, rd = UField.getLongitudePlaneB(5, 175*np.pi/180, 10)

# UField.nMax = true_nMax
# Bc, rc = UField.getLongitudePlaneB(5, 175*np.pi/180, 10)

# ratios = np.linalg.norm(Bc - Bd, axis = 2)/np.linalg.norm(Bd)
# ratios = np.nan_to_num(ratios)
# print(ratios)
# maxratio = np.amax(ratios, axis = 1)
# maxratio = np.amax(maxratio, axis = 0)
# print(maxratio)


# x = r[:, :, 0]
# y = r[:, :, 1]
# z = r[:, :, 2]
# u = B[:, :, 0]
# v = B[:, :, 1]
# w = B[:, :, 2]

# BMag = np.linalg.norm(B, axis=2)
# print(BMag)

# xH, yH, zH, uH, vH, wH = tools.loadBField("Output/complete_field_phi=137_CI.npz")

# print(uH[0, 0]**2 + vH[0, 0]**2 + wH[0, 0]**2)
# print(u[0, 0]**2 + v[0, 0]**2 + w[0, 0]**2)

# print(np.shape(uH))
# print(np.shape(u))
# print(u - uH)
# print(uH)



# x = np.array([[0, 5, 10], [0, 5, 10], [0, 5, 10]])
# y = np.array([[0, 0, 0], [5, 5, 5], [10, 10, 10]])
# x = np.array([0, 5, 10])
# y = np.array([0, 5, 10])
# x = x - 2.5
# y = y - 2.5
# c = np.array([[1, 2], [3, 4]])

# plt.figure(6)
# plt.pcolormesh(x, y, c,cmap='plasma')




# plt.show()









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

#========================================================================================
# UField.plotDeviationColourMapLShell(1.5)
# UField.plotDeviationColourMapLShell(2)
# UField.plotDeviationColourMapLShell(2.5)
# UField.plotDeviationColourMapLShell(3)
# UField.plotDeviationColourMapLShell(3.5)
# UField.plotDeviationColourMapLShell(4)

UField.plotDeviationColourMapLongitudePlane([0, 10, 20, 30], 5, 50)

plt.show()