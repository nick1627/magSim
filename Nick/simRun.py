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
This will be a test of the simulation over a range of interest to check
our simulation is valid.  Not working!!!
"""
# Ru = 25600000           #radius of Uranus in metres

# #Create field of Uranus
# g = np.array([[11278, 10928, 0], [-9648, -12284, 1453]]) #these are in nanoteslas
# h = np.array([[0, -16049, 0], [0, 6405, 4220]])
# g = g/1000000000
# h = h/1000000000

# UField = SHField(Ru, g, h, 0, 0)
# UField.rotate("Field") #want field aligned coordinates

# BMag = np.linalg.norm(UField.getField(np.array([4*Ru, 0, 0])))
# uniformB = UniformField(np.array([0, 0, BMag]))

# initialPosition = np.array([0, 0, 0])
# initialVelocityDirection = np.array([1, 1, 1])
# initialKE = 10000 #eV

# # e = Electron(initialPosition, initialVelocityDirection, initialKE)
# p1 = Proton(initialPosition, initialVelocityDirection, initialKE)
# p2 = Proton(initialPosition, initialVelocityDirection, 10*initialKE)
# p3 = Proton(initialPosition, initialVelocityDirection, 100*initialKE)
# p4 = Proton(initialPosition, initialVelocityDirection, 1000*initialKE)



# protonList = [p1, p2, p3, p4]

# manager = SimulationManager(uniformB, protonList, 0.02, N=4, fileKeyWord="varyEnergy", endTimeList=10000)
# manager.runAllSims()
# manager.plotAllEnergy()


#===============================================================================================

# Ru = 25600000           #radius of Uranus in metres

# #Create field of Uranus
# g = np.array([[11278, 10928, 0], [-9648, -12284, 1453]]) #these are in nanoteslas
# h = np.array([[0, -16049, 0], [0, 6405, 4220]])
# g = g/1000000000
# h = h/1000000000

# UFieldWhole = SHField(Ru, g, h, 0, 0)
# UFieldDipole = copy.deepcopy(UFieldWhole)
# UFieldWhole.rotate("Field") #want field aligned coordinates
# UFieldDipole.rotate("Field")

# UFieldDipole.setDipoleOnly(True)

# initialPosition = np.array([6*Ru, 0, 0])
# initialVelocityDirection = np.array([0.1, 0.1, 1])
# initialKE = 10000 #eV

# # e = Electron(initialPosition, initialVelocityDirection, initialKE)
# p1 = Proton(initialPosition, initialVelocityDirection, initialKE)

# FieldList = [UFieldDipole]

# manager = SimulationManager(FieldList, p1, 0.02, N=1, fileKeyWord="nick-Uranus-DipoleOnly", endStepList=100000)
# manager.runAllSims()

#===============================================================================================

# Ru = 25600000           #radius of Uranus in metres

# #Create field of Uranus
# g = np.array([[11278, 10928, 0], [-9648, -12284, 1453]]) #these are in nanoteslas
# h = np.array([[0, -16049, 0], [0, 6405, 4220]])
# g = g/1000000000
# h = h/1000000000

# UFieldWhole = SHField(Ru, g, h, 0, 0, False)
# UFieldDipole = SHField(Ru, g, h, 0, 0, True)
# UFieldWhole.rotate("Field") #want field aligned coordinates
# UFieldDipole.rotate("Field")


# # BMag = np.linalg.norm(UFieldDipole.getField(np.array([6*Ru, 0, 0])))
# # uniformB = UniformField(np.array([0, 0, BMag]))

# initialPosition = np.array([6*Ru, 0, 0])
# initialVelocityDirection1 = np.array([1, 1, 1])
# initialVelocityDirection2 = np.array([0.1, 0.1, 1])
# initialKE = 10000 #eV

# # e = Electron(initialPosition, initialVelocityDirection, initialKE)
# p1 = Proton(initialPosition, initialVelocityDirection1, initialKE)
# p2 = Proton(initialPosition, initialVelocityDirection2, initialKE)

# FieldList = [UFieldDipole]

# manager = SimulationManager(UFieldDipole, [p1], 50, N=1, fileKeyWord="nick-Uranus-2ndAttempt-DipoleOnly-1:1:1", endStepList=4000000)
# manager.runAllSims()
# manager.plotAllEnergy()
# manager.simulations[0].plotPositionOnTime(z=True)

# #===============================================================================================

# Ru = 25600000           #radius of Uranus in metres

# #Create field of Uranus
# g = np.array([[11278, 10928, 0], [-9648, -12284, 1453]]) #these are in nanoteslas
# h = np.array([[0, -16049, 0], [0, 6405, 4220]])
# g = g/1000000000
# h = h/1000000000

# UFieldDipole = SHField(Ru, g, h, 0, 0, True)
# UFieldDipole.rotate("Field")

# initialPosition = np.array([6*Ru, 0, 0])
# initialVelocityDirection = np.array([0.1, 0.1, 1])
# initialKE = 10000 #eV

# p1 = Proton(initialPosition, initialVelocityDirection, initialKE)

# sim1 = Simulation(UFieldDipole, p1)
# sim1.run(10)
# print(sim1.position)

# manager = SimulationManager(UFieldDipole, [p1], 50, N=1, fileKeyWord="nick-Uranus-FullField-1:1:1", endStepList=4000000)
# manager.runAllSims()
# manager.plotAllEnergy()
# manager.simulations[0].plotPositionOnTime(z=True)

#===============================================================================================


# Ru = 25600000           #radius of Uranus in metres

# #Create field of Uranus
# g = np.array([[11278, 10928, 0], [-9648, -12284, 1453]]) #these are in nanoteslas
# h = np.array([[0, -16049, 0], [0, 6405, 4220]])
# g = g/1000000000 #now in teslas
# h = h/1000000000

# UFieldWhole = SHField(Ru, g, h, 0, 0, False)
# UFieldDipole = SHField(Ru, g, h, 0, 0, True)
# UFieldWhole.rotate("Field") #want field aligned coordinates
# UFieldDipole.rotate("Field")


# # BMag = np.linalg.norm(UFieldDipole.getField(np.array([6*Ru, 0, 0])))
# # uniformB = UniformField(np.array([0, 0, BMag]))

# # initialPosition = np.array([-1.879385242, -0.6840402867, 0])*Ru
# initialPosition = np.array([-1.215537244, 6.893654271, 0])*Ru
# # initialVelocityDirection1 = np.array([0.39511149, 0, 1])
# initialVelocityDirection1 = np.array([0.09350383699, 0, 1])
# initialKE = 1000000 #eV

# # e = Electron(initialPosition, initialVelocityDirection, initialKE)
# p1 = Proton(initialPosition, initialVelocityDirection1, initialKE)
# # p2 = Proton(initialPosition, initialVelocityDirection2, initialKE)


# manager = SimulationManager(UFieldWhole, [p1], 50, N=1, fileKeyWord="nick-Uranus-RegionTest2-FullField-L7", endStepList=1000000)
# manager.runAllSims()

#===============================================================================================




UFieldWhole = UranusField(False)
UFieldDipole = UranusField(True)
UFieldWhole.rotate("Field") #want field aligned coordinates
UFieldDipole.rotate("Field")


energyList = [2*10**7, 3*10**7, 4*10**7, 5*10**7, 6*10**7, 7*10**7, 8*10**7, 9*10**7] #energies in eV

manager1 = LocationCheck(7, 30, 200, 8, energyList, "proton", UFieldWhole, endStepList=250000, fileNameAddition="-Uranus-7-30-200-")
manager1.runAllSims()

# plt.show()