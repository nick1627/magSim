"""
This is a demonstration of the project code.

You will need to have Python 3.8.5 if you want to run it.  Dependencies can be found in requirements.txt at the top level of the repo.

WHAT IS THIS PROJECT ABOUT?
Charged particles (usually protons and electrons) get trapped in planetary magnetic fields in regions called radiation belts.
The distribution of these particles is unusual at Neptune and Uranus, and we think it's because of their magnetic fields, which
differ from the fields at most other planets in the solar system because they have a large 'quadrupole moment'.

To prove this, we make a simulation of the magnetic fields at Uranus using data from Voyager II.
Then we inject particles into this field and observe their subsequent motion.
Over long time scales, we show that radial motion exists, which would not exist for planets that have a dipole-only magnetic field.

The code will now be demonstrated.

Nick/fields.py contains the various field classes
Nick/simulation.py contains the simulation class and simulationManager class
Nick/particles.py contains the particle classes.
"""

#============================Import stuff (unimportant)==============================
# from math import dist
import sys, os
sys.path.insert(0, os.getcwd())

#First some imports
import numpy as np
import matplotlib.pyplot as plt
# import scipy as sp
# from scipy.optimize import curve_fit

#Import stuff that is mine
from fields import *
from simulation import *
from particles import *


#Useful to have radius of Uranus
Uradius = 25600000

#Create the magnetic field of Uranus (see the fields file)
UField = UranusField(False)
#Rotate reference frame, since the axis of the dipole field is not aligned with the 
#geographical pole of Uranus, just like on Earth.
UField.rotate("Field")


#Now create a proton
initialPosition = np.array([4*Uradius, 0, 0]) #initial position on the equatorial plane but 4 radii out
initialVelocityDirection = np.array([0, 0, 1]) # point the particle northwards
initialKE = 10**8 #Initial kinetic energy in electronvolts 

p = Proton(initialPosition, initialVelocityDirection, initialKE) #create proton object p

#Now insert these into the simulation
sim1 = Simulation(UField, p)

#run the simulation
sim1.run(endStep = 10000) #run for 10k steps.  This takes me about 12 seconds

sim1.plotPositionOnTime() #plot a track of the position over time

plt.show() #show the graphs.  Graphs may overlap


