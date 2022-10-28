"""
Generate all data/figures for report
"""

from re import S
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


#=======================================field plots==============================================
#Field plots
# UField1 = UranusField(dipoleOnly=True)
# UField2 = UranusField(dipoleOnly=False)
# UField1.rotate("Field")
# UField2.rotate("Field")
# UField1.plot2DField("x", 0, 16, 25600000/4, 25600000)


# UField2.plotDeviationColourMapLShell(7)
# UField2.plotLShellBMagnitude(2)
# UField2.plotDeviationColourMapLShell2([2, 7], save="Output/Figures/L2-7BDeviationNick2.eps")

#=======================================drift plots==============================================
# UField1 = UranusField(dipoleOnly=True)
# UField1.rotate("Field")
# # UField1.plotDriftDirectionLShell([7])

UField2 = UranusField(dipoleOnly=False)
UField2.rotate("Field")
UField2.plotDriftDirectionLShell([2,7])#, save="Output/Figures/L2-7DriftNick2.eps")



#======================================3D position plot======================================
# sim3 = Simulation(simDataPath="Output/FinalData/Nick/MainGraph/locationCheck-Uranus-gammaRemoved-7-30-200-0--Proton-dipoleOnly-1000000.npz")
# sim3.plotPositionOnTime()



#======================================Constants plot======================================

# sim4 = Simulation(simDataPath="Output/locationCheckLong-constantPlot-50Step-7-40-0-0--Proton-dipoleOnly-10000000.npz")


# sim4.plotConstantsOnTime(save = "Output/Figures/constants3Nick.eps")


#=====================================Main result 1===============================================

# Uradius = 25600000
# data = tools.loadRegionData("Output/RegionTests/NICK_ONLY_regionTest_Uranus_7-30-200.npz")
# tools.plotRChangeOnEnergy2(data, Uradius, 7, 30, 200, save="Output/Figures/mainResultNick.eps")

# #=====================================Main result 2===============================================

# tools.plotDeltaRHistogram("Output/RegionTests/randomiseB_Nick_regionTest_Uranus_7-30-200.npz", 25600000, save="Output/Figures/randomiseBHistogramNick.eps")


# #=====================================Main result 3===============================================


# tools.plotCircumferenceGraphs("Output/RegionTests/NICK_ONLY_regionTest_CircumferenceNorth_Uranus_7-30-200.npz", "Output/RegionTests/NICK_ONLY_regionTest_CircumferenceSouth_Uranus_7-30-200.npz", 25600000, save="Output/Figures/circumferenceNick.eps")

# #=====================================final thing===============================================


# UField2 = UranusField(dipoleOnly=False)
# UField2.rotate("Field")
# UField2.plotDriftDirectionLShell([7], particlePath = "Output/FinalData/Nick/CircumferenceNorth/locationCheck-Uranus-Circumference-N220-50Step-7-30-220-0--Proton-fullField-1000000.np")#, save="Output/Figures/L2-7DriftNick2.eps")


plt.show()