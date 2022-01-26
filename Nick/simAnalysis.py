"""
Here we analyse the results of the simulations
"""

#First some imports
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import curve_fit

#Import stuff that is mine
from fields import *
from simulation import *
from particles import *
#=======================================================================================
#Analyse results of full field vs dipole only field
dipoleSim = Simulation(simDataPath="Output/nick-Uranus-DipoleOnly-Proton-10000.npz")
fullSim = Simulation(simDataPath="Output/nick-Uranus-FullField-Proton-10000.npz")


# dipoleSim.plotPositionOnTime()
# fullSim.plotPositionOnTime()

dipoleSim.plotLShellOnTime("- Dipole Only")
fullSim.plotLShellOnTime("- Full Field")

# dipoleSim.plotKEOnTime("- Dipole Only")
fullSim.plotKEOnTime("- Full Field")

# fullSim.plotFirstAIOnTime()
dipoleSim.plotFirstAIOnTime()

plt.show()