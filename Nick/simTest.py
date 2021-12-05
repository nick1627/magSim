"""
This is the testing script for the simulation
"""
#modify path so we get access to the stuff we need
import sys, os
sys.path.insert(0, os.getcwd())

#Import stuff that isn't mine
import numpy as np
import matplotlib.pyplot as plt

#Import stuff that is mine
from fields import *
from simulation import *
from particles import *
import tools


#===============================================================================
uniformB = UniformField(np.array([0, 0, 1]))
initialPosition = np.array([0, 0, 0])
initialVelocity = np.array([0.01, 0.01, 0.01])
e = Electron(initialPosition, initialVelocity)
sim = Simulation(uniformB, e, 0.02)
sim.run(50, returnData=False, saveData=False, naturalUnits=True)
sim.plotPositionOnTime()
sim.plotKEOnTime()



#===============================================================================
plt.show()
