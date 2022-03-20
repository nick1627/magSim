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
import dataAnalysis
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

# Uradius = 25600000

# sim = Simulation(simDataPath = "Output/Trash/locationCheck-Uranus-gammaRemoved-7-30-200-0--Proton-fullField-10000000.npz")
# i = sim.getFirstEquatorialIndex()
# v = sim.velocity[i]
# deltaT = sim.time[i] - sim.time[i-1]
# print(np.linalg.norm(v*deltaT/Uradius))
# # sim.plotAltitudeOnTime()
# # sim.plotPositionOnTime(z=True)
# sim.plotKEOnTime()
# dataAnalysis.findFirstEnergyStepHeight(sim)
# # sim.plotFirstAIOnTime()
# # sim1.plotLShellOnTime()
# # sim.saveBounceData("Output/RegionTests/regionTest_Uranus_7-30-200.npz")
# # print(tools.loadRegionData("Output/RegionTests/regionTest_Uranus_7-30-200.npz"))
# plt.show()


# Uradius = 25600000
# data = tools.loadRegionData("Output/RegionTests/regionTest_Uranus_7-30-200.npz")
# fullData = tools.selectCriteria(data, species="proton", field="fullField")
# dipoleData = tools.selectCriteria(data, species="proton", field="dipoleOnly")
# tools.plotRChangeOnEnergy(dipoleData, Uradius, 7, 30, 200)
# tools.plotRChangeOnEnergy(fullData, Uradius, 7, 30, 200)
# plt.show()

# tools.deleteOlderThan(218, "Nick", "Output/RegionTests/regionTest_Uranus_7-30-200.npz")

# Uradius = 25600000
# data = tools.loadRegionData("Output/RegionTests/regionTest_Uranus_7-30-200.npz")
# fullData = tools.selectCriteria(data, species="proton", field="fullField")
# dipoleData = tools.selectCriteria(data, species="proton", field="dipoleOnly")
# tools.plotRChangeOnEnergy(fullData, Uradius, 7, 30, 200)
# tools.plotRChangeOnEnergy(dipoleData, Uradius, 7, 30, 200)
# tools.plotGyroradiusOnEnergy(fullData)


# Uradius = 25600000
# data = tools.loadRegionData("Output/RegionTests/regionTest_Uranus_7-30-200.npz")
# tools.plotRChangeOnEnergy2(data, Uradius, 7, 30, 200)
# plt.show()


# long run plot ========================================================================
# sim1 = Simulation(simDataPath = "Output/Trash/locationCheck-Uranus-LongRunPlot-7-30-200-0--Proton-dipoleOnly-100000.npz")
# sim2 = Simulation(simDataPath = "Output/Trash/locationCheck-Uranus-LongRunPlot-7-30-200-0--Proton-fullField-100000.npz")

# sim1.plotLShellOnTime(otherSims = [sim2], legendList = ["Dipole only", "Complete field"], titleAddition=" - Targeting L=7, " + r"$\theta$" + "=30, " + r"$\phi$" + "=200, 0 phase")
# plt.savefig("Output/Figures/longRunLShell.eps", type="eps")
# plt.show()

# larmor error ========================================================================

#dipole only, 0 phase runs

# sim = Simulation(simDataPath = "Output/locationCheck-Uranus-50Step-7-30-200-0--Proton-fullField-10000.npz")
# sim.plotPositionOnTime(z=True)
# equatorialIndex = sim.getFirstEquatorialIndex()
# print(equatorialIndex)
# # deltaLarmorRadius = sim.getChangeInLarmorRadius(equatorialIndex, sim.particle.initialEnergy)
# # print("change in larmor radius/a")
# # print(deltaLarmorRadius/Uradius)
# # print("final larmor radius")
# # print(sim.getLarmorRadius(equatorialIndex)/Uradius)
# # print("Error in calculating guiding centre position")
# # print(sim.getEquatorialGuidingCentreError()/Uradius)
# # sim.plotAltitudeOnTime()
# # sim.plotPositionOnTime(z=True)
# # sim.plotKEOnTime()
# # sim.plotFirstAIOnTime()
# # sim1.plotLShellOnTime()
# # sim.saveBounceData("Output/RegionTests/regionTest_Uranus_7-30-200.npz")
# # print(tools.loadRegionData("Output/RegionTests/regionTest_Uranus_7-30-200.npz"))
# error = sim.getPositionAndVelocityError(equatorialIndex)
# # print(error)
# print(error[0]/Uradius)
# plt.show()



#======================analysing data for 800 steps=================================================


# sim = Simulation(simDataPath = "Output/locationCheck-Uranus-800Step-7-30-200-0--Proton-fullField-10000000.npz")
# sim.saveBounceData("Output/RegionTests/regionTest_Uranus_MoreSteps_7-30-200.npz")

# Uradius = 25600000
# data = tools.loadRegionData("Output/RegionTests/regionTest_Uranus_MoreSteps_7-30-200.npz")
# tools.plotRChangeOnEnergy2(data, Uradius, 7, 30, 200)
# plt.show()



#======================================Collecting Nick only results==========================================
# dataAnalysis.collectResults("Output/FinalData/Nick", "Output/RegionTests/NICK_ONLY_regionTest_Uranus_7-30-200.npz")

# Uradius = 25600000
# data = tools.loadRegionData("Output/RegionTests/NICK_ONLY_regionTest_Uranus_7-30-200.npz")
# tools.plotRChangeOnEnergy2(data, Uradius, 7, 30, 200)
# plt.show()


#The following files are ones that don't cross the equator

# sim = Simulation(simDataPath="Output/locationCheck-Uranus-200Step-7-30-200-90--Electron-dipoleOnly-10000000.npz")
# sim.plotPositionOnTime(z=True)
# plt.show()
# "Output/FinalData/Nick/locationCheck-Uranus-50Step-7-30-200-90--Proton-fullField-1000.npz" done checked
# "Output/FinalData/Nick/locationCheck-Uranus-7-30-200-270--Electron-fullField-100000.npz" done checked BUT MAY HAVE TO REDO AS ONLY USED 50 STEPS INSTEAD OF LIKE 200
# "Output/FinalData/Nick/locationCheck-Uranus-50Step-7-30-200-270--Proton-dipoleOnly-1000.npz" done
# "Output/FinalData/Nick/locationCheck-Uranus-gammaRemoved-7-30-200-90--Proton-fullField-4000000.npz" done checked
# "Output/FinalData/Nick/locationCheck-Uranus-50Step-7-30-200-180--Proton-dipoleOnly-1000.npz" done
# "Output/FinalData/Nick/locationCheck-Uranus-50Step-7-30-200-0--Proton-dipoleOnly-1000.npz" done checked
# "Output/FinalData/Nick/locationCheck-Uranus-50Step-7-30-200-180--Proton-fullField-1000.npz" done
# "Output/FinalData/Nick/locationCheck-Uranus-50Step-7-30-200-90--Proton-dipoleOnly-1000.npz" done
# "Output/FinalData/Nick/locationCheck-Uranus-gammaRemoved-7-30-200-0--Proton-dipoleOnly-6000000.npz" done


#also missing 10**7 electrons

#======================================checking errors==========================================

# sim = Simulation(simDataPath="Output/FinalData/Nick/locationCheck-Uranus-50Step-7-30-200-0--Proton-fullField-1000.npz")
# equatorialIndex = sim.getFirstEquatorialIndex()
# positionError, velocityError = sim.getPositionAndVelocityError(equatorialIndex)
# print(positionError/Uradius)

#======================================general analysis ==========================================

# sim = Simulation(simDataPath="Output/locationCheck-Uranus-Circumference-S56-50Step-7-150-56-0--Proton-fullField-1000000.npz")
# sim.plotPositionOnTime(z=True)
# plt.show()

# dataAnalysis.collectResults("Output/FinalData/Nick/CircumferenceSouth", "Output/RegionTests/NICK_ONLY_regionTest_CircumferenceSouth_Uranus_7-30-200.npz")

# tools.plotCircumferenceGraphs("Output/RegionTests/NICK_ONLY_regionTest_CircumferenceNorth_Uranus_7-30-200.npz", "Output/RegionTests/NICK_ONLY_regionTest_CircumferenceSouth_Uranus_7-30-200.npz", 25600000)

# plt.show()


#=====================================nice 3d plot==========================================

# sim = Simulation(simDataPath="Output/locationCheckdemoPic-50Step-5-50-0-0--Proton-dipoleOnly-10000000.npz")
# sim.plotPositionOnTime()
# plt.show()