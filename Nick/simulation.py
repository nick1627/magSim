"""
This file stores the code that analyses the simulation
"""
import numpy as np
from fields import *
from particles import *

class Simulation:
    def __init__(self, field, particle, timestep):
        self.timeStep = timestep
        self.field = field
        self.particle = particle

        #Diagnostic info
        self.positions = []
        self.velocities = []

    def run(self, endTime, returnData=True, saveData=True):
        currentTime = 0
        #Record initial values
        self.positions.append(self.particle.getPosition())
        self.velocities.append(self.particle.getVelocity())

        while currentTime < endTime:
            #update particle position
            self.particle.updatePosition(self.field, self.timeStep)
            #Record data
            self.positions.append(self.particle.getPosition())
            self.velocities.append(self.particle.getVelocity())
            #Now increment time for the next step
            currentTime += self.timeStep
        
        #Simulation has now completed
        if saveData:
            self.positions = np.array(self.positions)
            self.velocities = np.array(self.velocities)

        if returnData:
            return self.positions, self.velocities
        else:
            return
    
    def plotLShellOnTime(self):
        return
  
