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
sim = Simulation(1)
