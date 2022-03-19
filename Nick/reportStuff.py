"""
Generate all data/figures for report
"""

import sys, os
sys.path.insert(0, os.getcwd())

#First some imports
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import curve_fit

#Import stuff that is mine
from fields import *
from simulation import *
from particles import *
import tools
import dataAnalysis


#=======================================drift plots==============================================
# UField1 = UranusField(dipoleOnly=True)
# UField1.rotate("Field")
# # UField1.plotDriftDirectionLShell([7])

UField2 = UranusField(dipoleOnly=False)
UField2.rotate("Field")
UField2.plotDriftDirectionLShell([2,7], save="Output/Figures/L2-7DriftNick.png")




# # sim1 = Simulation(simDataPath="Output/FinalData/Nick/MainGraph/locationCheck-Uranus-gammaRemoved-7-30-200-0--Proton-dipoleOnly-1000000.npz")
# sim2 = Simulation(simDataPath="Output/FinalData/Nick/CircumferenceSouth/locationCheck-Uranus-Circumference-S282-50Step-7-150-282-0--Proton-fullField-1000000.npz")
# index = sim2.getFirstEquatorialIndex()
# sim2.plotLShellOnTime(finalIndex=index)

# # # firstPosition = sim.position[0, :]
# # # print(firstPosition)

# # # sim.plotPositionOnTime()
# # # vec = sim.convertCartesianToPolar(firstPosition)
# # # vec[1]*=180/np.pi
# # # vec[2]*=180/np.pi
# # # vec[2] = vec[2] % 360
# # # print(vec)

# # path1 = sim1.getFirstBounceThetaPhi()
# path2 = sim2.getFirstBounceThetaPhi()

# # # UField1.plotDriftDirectionLShell([7], particlePath = path1)
# UField2.plotDriftDirectionLShell([7], particlePath = path2)
# # UField2.plotDriftDirectionLongitudePlane(phi=200, rMax=3, N=200)

# # plt.show()


# def normVec(r):
#     x = r[0]
#     y = r[1]
#     z = r[2]

#     xc = x - 2*(x*z**2)/(x**2 + y**2)
#     yc = y - 2*(y*z**2)/(x**2 + y**2)
#     zc = 3*z

#     vec = np.array([xc, yc, zc])
#     return vec/np.linalg.norm(vec)

# theta = np.linspace(0, np.pi, 180)
# r = 2*np.sin(theta)**2

# L = np.zeros((np.shape(r)[0], 3))
# for i in range(0, len(theta)):
#     L[i, 2] = r[i]*np.cos(theta[i])
#     L[i, 0] = r[i]*np.sin(theta[i])

# norms = np.zeros(np.shape(L))
# for i in range(0, len(theta)):
#     norms[i, :] = normVec(L[i, :])






# plt.figure(1)
# plt.plot(L[:,0], L[:,2])
# for i in range(0, len(theta)):
#     plt.arrow(L[i,0], L[i,2], norms[i,0], norms[i,2], head_width=0.02)
# plt.show()

#work in y=0 plane

# print(UField1.getGradB(np.array([200000, 0, 0])))


#=====================================Main result 1===============================================

# Uradius = 25600000
# data = tools.loadRegionData("Output/RegionTests/NICK_ONLY_regionTest_Uranus_7-30-200.npz")
# tools.plotRChangeOnEnergy2(data, Uradius, 7, 30, 200, save="Output/Figures/mainResultNick.eps")

# #=====================================Main result 2===============================================

# tools.plotDeltaRHistogram("Output/RegionTests/randomiseB_Nick_regionTest_Uranus_7-30-200.npz", 25600000, save="Output/Figures/randomiseBHistogramNick.eps")


# #=====================================Main result 3===============================================


# tools.plotCircumferenceGraphs("Output/RegionTests/NICK_ONLY_regionTest_CircumferenceNorth_Uranus_7-30-200.npz", "Output/RegionTests/NICK_ONLY_regionTest_CircumferenceSouth_Uranus_7-30-200.npz", 25600000, save="Output/Figures/circumferenceNick.eps")



plt.show()