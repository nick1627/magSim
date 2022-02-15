"""
Here we analyse the results of the simulations
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
#=======================================================================================
#Analyse results of full field vs dipole only field
#dipoleSim = Simulation(simDataPath="Output/nick-Uranus-DipoleOnly-Proton-10000.npz")
# dipoleSim = Simulation(simDataPath="Output/Trash/nick-Uranus-RegionTest1-DipoleOnly--Proton-100000.npz")
# fullSim = Simulation(simDataPath="Output/Trash/nick-Uranus-RegionTest2-FullField-L7-Proton-1000000.npz")

# # dipoleSim.plotPositionOnTime()
# fullSim.plotPositionOnTime()

# # dipoleSim.plotLShellOnTime("- Dipole Only")
# fullSim.plotLShellOnTime("- Full Field")

# # dipoleSim.plotKEOnTime("- Dipole Only")
# fullSim.plotKEOnTime("- Full Field")

# # dipoleSim.plotFirstAIOnTime("- Dipole Only")
# fullSim.plotFirstAIOnTime("- Full Field")

# print("plotting now")
# plt.show()
# print("plotted")


#=======================================================================================


# sim = Simulation(simDataPath = "Output/Trash/locationCheck-Uranus-7-30-200-90--Proton-fullField-1000000.npz")
# sim.plotPositionOnTime(z=True)
# sim.plotKEOnTime()
# sim.plotFirstAIOnTime()
# sim.plotLShellOnTime()
# sim.saveBounceData("Output/RegionTests/regionTest_Uranus_7-30-200.npz")
# print(tools.loadRegionData("Output/RegionTests/regionTest_Uranus_7-30-200.npz"))

# Uradius = 25600000
# data = tools.loadRegionData("Output/RegionTests/regionTest_Uranus_7-30-200.npz")
# fullData = tools.selectCriteria(data, species="proton", field="fullField")
# dipoleData = tools.selectCriteria(data, species="proton", field="dipoleOnly")
# tools.plotRChangeOnEnergy(fullData, Uradius, 7, 30, 200)
# tools.plotRChangeOnEnergy(dipoleData, Uradius, 7, 30, 200)

Uradius = 25600000
data = tools.loadRegionData("Output/RegionTests/regionTest_Uranus_7-30-200.npz")
fullData = tools.selectCriteria(data, species="proton", field="fullField")
# dipoleData = tools.selectCriteria(data, species="proton", field="dipoleOnly")
# tools.plotRChangeOnEnergy(fullData, Uradius, 7, 30, 200)
# tools.plotRChangeOnEnergy(dipoleData, Uradius, 7, 30, 200)
tools.plotGyroradiusOnEnergy(fullData)


plt.show()