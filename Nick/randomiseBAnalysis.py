"""
In this script we randomise the quadrupole component of the B field within the range allowed by the errors.
"""

from math import dist
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



# g = np.array([[11278, 10928, 0], [-9648, -12284, 1453]]) #these are in nanoteslas
# h = np.array([[0, -16049, 0], [0, 6405, 4220]])

# g_error_q = np.array([[0, 0, 0], [352, 425, 709]])
# h_error_q = np.array([[0, 0, 0], [0, 541, 578]])

# g = g/1000000000 #now in teslas
# h = h/1000000000

# g_error_q = g_error_q/1000000000
# h_error_q = h_error_q/1000000000

# g_origin = g - g_error_q
# h_origin = h - h_error_q

# delta_r = []

# print(np.shape(tools.loadRegionData("Output/RegionTests/randomiseB_regionTest_Uranus_7-30-200.npz"))[0])

# for i in range(0, 694):
#     #First create the field objects
#     randomNumbers_g = np.random.rand(np.shape(g)[0], np.shape(g)[1])
#     randomNumbers_h = np.random.rand(np.shape(h)[0], np.shape(h)[1])
    
#     current_g = g_origin + 2*randomNumbers_g*g_error_q
#     current_h = h_origin + 2*randomNumbers_h*h_error_q

#     UField = SHField(25600000, current_g, current_h, dipoleOnly=False)
#     UField.rotate("Field")

#     #Now do the runs
#     #Don't forget to reset it to 50 step!!!
#     manager1 = LocationCheck(7, 30, 200, 1, [10**7], "proton", 0, UField, endStepList=16000, fileNameAddition="-Uranus-randomB-50Step-")
#     manager1.runAllSims()

#     #Remove these lines when ready
#     # manager1.simulations[0].plotPositionOnTime(z=True)
#     # plt.show()
#     try:
#         manager1.simulations[0].saveBounceData("Output/RegionTests/randomiseB_regionTest_Uranus_7-30-200.npz")
#     except:
#         manager1.simulations[0].plotPositionOnTime(z=True)


tools.plotDeltaRHistogram("Output/RegionTests/randomiseB_Nick_regionTest_Uranus_7-30-200.npz", 25600000)
plt.show()


