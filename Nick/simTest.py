"""
This is the testing script for the simulation
"""
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
import tools


#===============================================================================

# Ru = 25600000           #radius of Uranus in metres

# g = np.array([[11278, 10928, 0], [-9648, -12284, 1453]]) #these are in nanoteslas
# h = np.array([[0, -16049, 0], [0, 6405, 4220]])
# g = g/1000000000
# h = h/1000000000

# UField = SHField(Ru, g, h, 0, 0)
# UField.rotate("Field")

# BMag = np.linalg.norm(UField.getField(np.array([6*Ru, 0, 0])))
# BMag = 1
# uniformB = UniformField(np.array([0, 0, BMag]))

# initialPosition = np.array([0, 0, 0])
# initialVelocityDirection = np.array([1, 1, 1])
# initialKE = 8.54E-9#100 #keV

# # e = Electron(initialPosition, initialVelocityDirection, initialKE)
# p = Proton(initialPosition, initialVelocityDirection, initialKE)

# sim = Simulation(uniformB, p, 0.02)
# sim.run(endOnTime=False, endTime = 1000, endStep=10000, naturalUnits=True)
# sim.plotPositionOnTime()
# sim.plotKEOnTime()
# sim.plotVelocityOnTime()
# sim.plotVelocityErrorOnTime()

# plt.figure(10)
# startEnergy = np.array([1, 5, 10, 20, 30, 100, 200])
# percentageGain = np.array([-11.6, -11.6, -11.6, -11.5, -11.4, -10.9, -10.3])
# plt.plot(startEnergy, percentageGain, marker="x", color="black", linestyle="none")
# plt.title("Percentage gain in energy against initial energy")
# plt.xlabel("Initial energy (keV)")
# plt.ylabel("Percentage gain in energy over simulation")

# plt.figure(20)
# timeStep = np.array([0.001, 0.005, 0.010, 0.015, 0.020, 0.030, 0.040])
# percentageGain = np.array([-0.427, -2.119, -4.198, -6.235, -8.23, -12.106, -15.834])
# plt.plot(timeStep, percentageGain, marker="x", color="black", linestyle="none")
# plt.title("Percentage gain in energy against time step for a constant total run time (50 gyrations ish)")
# plt.xlabel("Time steps (fractions of initial gyro-period)")
# plt.ylabel("Percentage gain in energy over simulation")

# plt.figure(30)
# timeStep = np.array([0.001, 0.005, 0.010, 0.015, 0.020, 0.030, 0.040])
# percentageGain = np.array([-0.009, -0.214, -0.853, -1.909, -3.371, -7.435, -12.857])

# x = np.linspace(0, timeStep[-1], 1000)
# def func(x, A, B):
#     return B*x**A

# params, cov = curve_fit(func, timeStep, percentageGain, bounds = ([1, -10000], [10, -1]))

# plt.plot(timeStep, percentageGain, marker="x", color="black", linestyle="none")
# plt.plot(x, func(x, params[0], params[1]), label="y = " + str(params[1]) + "x^" + str(params[0]))
# plt.title("Percentage gain in energy against time step for 1000 steps")
# plt.xlabel("Time step (fractions of initial gyro-period)")
# plt.ylabel("Percentage gain in energy over simulation")
# plt.legend()

# plt.figure(40)
# timeStep = np.array([0.001, 0.005, 0.010, 0.015, 0.020, 0.030, 0.040, 0.05, 0.07, 0.1])#, 0.2])
# percentageGain = np.array([-0.001, -0.018, -0.074, -0.167, -0.299, -0.69, -1.269, -2.069, -4.518, -11.191])#, -40.68])

# x = np.linspace(0, timeStep[-1], 1000)
# def func(x, B, C):
#     return B*x**6 + C*x**0.5

# def func2(x, B, C, D):
#     return B*C**x + D


# params, cov = curve_fit(func2, timeStep, percentageGain, bounds = ([-10000, 1, 0], [0, np.inf, 4]))

# plt.plot(timeStep, percentageGain, marker="x", color="black", linestyle="none")
# # plt.plot(x, func2(x, params[0], params[1], params[2]), label="y = " + str(params[1]) + "x^" + str(params[0]))
# plt.title("Percentage gain in speed against time step for 1000 steps")
# plt.xlabel("Time step (fractions of initial gyro-period)")
# plt.ylabel("Percentage gain in speed over simulation")
# plt.legend()


#==================start testing of simulation manager==========================

# Ru = 25600000           #radius of Uranus in metres

# g = np.array([[11278, 10928, 0], [-9648, -12284, 1453]]) #these are in nanoteslas
# h = np.array([[0, -16049, 0], [0, 6405, 4220]])
# g = g/1000000000
# h = h/1000000000

# UField = SHField(Ru, g, h, 0, 0)
# UField.rotate("Field")

# BMag = np.linalg.norm(UField.getField(np.array([6*Ru, 0, 0])))
# #BMag = 1
# uniformB = UniformField(np.array([0, 0, BMag]))

# initialPosition = np.array([0, 0, 0])
# initialVelocityDirection = np.array([1, 1, 1])
# initialKE = 10000 #eV

# # e = Electron(initialPosition, initialVelocityDirection, initialKE)
# p1 = Proton(initialPosition, initialVelocityDirection, initialKE)
# # p2 = Proton(initialPosition, initialVelocityDirection, 2*initialKE)
# # p3 = Proton(initialPosition, initialVelocityDirection, 3*initialKE)
# # p4 = Proton(initialPosition, initialVelocityDirection, 4*initialKE)
# # p5 = Proton(initialPosition, initialVelocityDirection, 5*initialKE)
# protonList = [p1]

# manager = SimulationManager(uniformB, protonList, 50, N=1, fileKeyWord="newSimTest", endStepList=1000)
# manager.runAllSims()






#===============================================================================
# Ru = 25600000           #radius of Uranus in metres



# UField = UranusField(False)
# UField.rotate("Field")

# # print(UField.getField([6*Ru, 0, 0]))

# BMag = np.linalg.norm(UField.getField(np.array([6*Ru, 0, 0])))
# # #BMag = 1
# uniformB = UniformField(np.array([0, 0, BMag]))

# initialPosition = np.array([0, 0, 0])
# initialVelocityDirection = np.array([1, 0, 0])
# initialKE = 10**9 #eV

# p1 = Proton(initialPosition, initialVelocityDirection, initialKE)

# s1 = Simulation(uniformB, p1)
# s1.run(50000)

# # s1.plotPositionOnTime()
# larmorchange = s1.getLarmorRadius(-1) - s1.getLarmorRadius(0)
# print(larmorchange)
#===============================================================================
#===============================================================================
Ru = 25600000           #radius of Uranus in metres



UField = UranusField(False)
UField.rotate("Field")

# print(UField.getField([6*Ru, 0, 0]))

# BMag = np.linalg.norm(UField.getField(np.array([6*Ru, 0, 0])))
# #BMag = 1
# uniformB = UniformField(np.array([0, 0, BMag]))

# initialPosition = np.array([-168458653.14769566, -61313935.450335704, 1.0972835320360285e-08])
# initialVelocityDirection = np.array([1382649.4412019397, -3798798.1187238824, 43234552.18773658])
initialPosition = np.array([-169052073.9254811, -61529922.94984471, 1.0972835320360285e-08])

initialVelocityDirection = np.array([1382649.4412019397, -3798798.1187238824, 43234552.18773658])

initialKE = 10**7 #eV

p1 = Proton(initialPosition, initialVelocityDirection, initialKE)

s1 = Simulation(UField, p1)
s1.run(endStep=7200)
print("sim complete")
lastPosition = s1.position[7200]

positions = s1.position
for i in range(0, np.shape(positions)[0]):
    if positions[i, 2] < 0 and i > 7000:
        flipindex = i
        break

print(flipindex)
thingy = positions[flipindex]
print(positions[flipindex])

print(lastPosition)
print("hi")

rL = s1.getLarmorRadius(flipindex)

difference = np.array([-1207972.22539455, -590354.15896964, 418455.22015373])
mag = np.linalg.norm(difference)
print(mag/rL)
print("larmor radius")
print(rL)

plt.show()


