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



#=====================================Main result 1===============================================

Uradius = 25600000
data = tools.loadRegionData("Output/RegionTests/NICK_ONLY_regionTest_Uranus_7-30-200.npz")
tools.plotRChangeOnEnergy2(data, Uradius, 7, 30, 200, save="Output/Figures/mainResultNick.eps")

#=====================================Main result 2===============================================

tools.plotDeltaRHistogram("Output/RegionTests/randomiseB_Nick_regionTest_Uranus_7-30-200.npz", 25600000, save="Output/Figures/randomiseBHistogramNick.eps")


#=====================================Main result 3===============================================


tools.plotCircumferenceGraphs("Output/RegionTests/NICK_ONLY_regionTest_CircumferenceNorth_Uranus_7-30-200.npz", "Output/RegionTests/NICK_ONLY_regionTest_CircumferenceSouth_Uranus_7-30-200.npz", 25600000, save="Output/Figures/circumferenceNick.eps")



plt.show()