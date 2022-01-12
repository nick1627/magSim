"""
This is the file where major non-test runs of the simulation go.

Difficult tests will also go here, like tests confirming the final range of the simulation
validity in terms of energy etc.

"""
#================
#Imports
#modify path so we get access to the stuff we need
import sys, os
sys.path.insert(0, os.getcwd())

#Import stuff that isn't mine
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import curve_fit

#Import stuff that is mine
from fields import *
from simulation import *
from particles import *

#======================================================================
"""
This will be a test of the simulation over a range of interest to check th
"""
Ru = 25600000           #radius of Uranus in metres

#Create field of Uranus
g = np.array([[11278, 10928, 0], [-9648, -12284, 1453]]) #these are in nanoteslas
h = np.array([[0, -16049, 0], [0, 6405, 4220]])
g = g/1000000000
h = h/1000000000

UField = SHField(Ru, g, h, 0, 0)
UField.rotate("Field") #want field aligned coordinates

BMag = np.linalg.norm(UField.getField(np.array([2*Ru, 0, 0])))
uniformB = UniformField(np.array([0, 0, BMag]))

initialPosition = np.array([0, 0, 0])
initialVelocityDirection = np.array([1, 1, 1])
initialKE = 10000000 #eV

# e = Electron(initialPosition, initialVelocityDirection, initialKE)
p1 = Proton(initialPosition, initialVelocityDirection, initialKE)

protonList = [p1]

manager = SimulationManager(uniformB, protonList, 0.000001, N=1, fileKeyWord="Major1", endStepList=1000000)
manager.runAllSims()
manager.plotAllEnergy()







plt.show()