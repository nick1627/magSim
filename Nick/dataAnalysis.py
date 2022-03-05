"""
Miscellaneous data analysis code goes here
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


def findFirstEnergyStepHeight(sim):
    """
    sim:        The simulation object with data loaded
    """

    E = sim.plotKEOnTime(plot=False, returnData=True)

    #first attempt uses hope that energy is initially constant
    # if E[1] != E[0]:
    #     raise(Exception("this method won't work"))
    # else:
    #     flat=True
    #     counter=0
    #     while flat==True:
    #         if E[counter+1] == E[counter]:
    #             counter += 1
    #         else:
    #             flat=False
        
    #     while flat==False:
    #         if E[counter+1] == E[counter]:
    #             flat = True
    #         else:
    #             counter += 1
        
    #     #should have reached top of step
    #     print(sim.time[counter])

    #Need to differentiate
    dEdt = []
    counter = 0
    for i in range(0, np.shape(E)[0]-1):
        dEdt.append((E[i+1] - E[i])/(sim.time[i+1] - sim.time[i]))

    dEdt = np.array(dEdt)

    ft = np.fft.rfft(dEdt)
    ft[0] = 0
    ft = abs(ft)
    maxFreqIndex = np.argmax(ft)
    
    freq = np.fft.rfftfreq(ft.size)
    f = freq[maxFreqIndex]
    T = 1/f
    print(T)
    plt.figure(1)
    plt.plot(ft)
    
    plt.figure(2)
    plt.plot(sim.time, E)
    plt.show()






# sim = Simulation(simDataPath = "Output/locationCheck-ErrorTest-Uranus-7-30-200-0--Proton-fullField-10000000.npz")
# findFirstEnergyStepHeight(sim)



