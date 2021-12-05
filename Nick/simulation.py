"""
This file stores the code that analyses the simulation
"""
import numpy as np
from fields import *
from particles import *

class Simulation:
    def __init__(self, field, particle, timestep):
        #timestep always input as fraction of initial period of gyroradius
        self.timeStep = timestep
        self.field = field
        self.particle = particle

        #Diagnostic info
        self.position = []
        self.velocity = []
        self.time = []

    def run(self, endTime, returnData=True, saveData=True, naturalUnits=True):
        #End time is always in units of time steps, not seconds
        #Find values for conversion factors -- needed for calculating duration of timeStep
        initialB = self.field.getField(self.particle.getPosition())
        self.particle.computeConversionFactors(initialB)

        if naturalUnits:
            #EndTime is in units of time steps, not seconds
            currentTime = 0
     
                
            #Set up the field properly for calculations
            if self.field.getUnitState() == False:
                self.field.setNaturalUnits(True, self.particle.q, self.particle.m0, self.particle.larmorPeriod)
            
            #Convert particle to use natural units
            self.particle.convertToNatural()

            #Record initial values 
            self.position.append(self.particle.getPosition())
            self.velocity.append(self.particle.getVelocity())
            self.time.append(currentTime)

            while currentTime < endTime:
                #increment time 
                currentTime += self.timeStep
                #update particle position
                self.particle.updatePositionN(self.field, self.timeStep)
                #Record data
                self.position.append(self.particle.getPosition())
                self.velocity.append(self.particle.getVelocity())
                self.time.append(currentTime)
                
            
            #Simulation has now completed
            if self.field.getUnitState == True:
                self.field.naturalUnits(False, self.particle.q, self.particle.m0, self.particle.larmorPeriod)

            self.position = np.array(self.position)
            self.velocity = np.array(self.velocity)
            self.time = np.array(self.time)
            
            #Return units to SI state
            self.velocity = self.velocity*self.particle.c
            self.position = self.position*self.particle.larmorRadius
            self.time = self.time*self.particle.larmorPeriod
        
        else: #Not using natural units
             #EndTime is in units of time steps, not seconds
            currentTime = 0

            #Record initial values 
            self.position.append(self.particle.getPosition())
            self.velocity.append(self.particle.getVelocity())
            self.time.append(currentTime)

            endTime *= self.particle.larmorPeriod

            while currentTime < endTime:
                #increment time 
                currentTime += self.timeStep*self.particle.larmorPeriod
                #update particle position
                self.particle.updatePositionSI(self.field, self.timeStep*self.particle.larmorPeriod)
                #Record data
                self.position.append(self.particle.getPosition())
                self.velocity.append(self.particle.getVelocity())
                self.time.append(currentTime)
                print(currentTime)
            
            #Simulation has now completed
            self.position = np.array(self.position)
            self.velocity = np.array(self.velocity)

     

        if saveData:
            raise(Exception("you havent added saving data yet"))

        if returnData:
            return self.position, self.velocity
        else:
            return
        

    
    def plotLShellOnTime(self):
        return

    def plotPositionOnTime(self):
        ax = plt.figure().add_subplot(projection = "3d")
        ax.plot(self.position[:,0], self.position[:,1], self.position[:,2])
        ax.set_title("Position of particle in 3D")

        ax1 = plt.figure().add_subplot()
        ax1.plot(self.time, self.position[:, 0])
        ax1.set_title("x on time")
        ax2 = plt.figure().add_subplot()
        ax2.plot(self.time, self.position[:, 1])
        ax2.set_title("y on time")
        ax3 = plt.figure().add_subplot()
        ax3.plot(self.time, self.position[:, 2])
        ax3.set_title("z on time")

        ax4 = plt.figure().add_subplot()
        ax4.plot(self.position[:, 0], self.position[:, 1])
        ax4.set_aspect("equal")
        ax4.set_title("x and y positions as time evolves")

    def plotKEOnTime(self):
        print("warning:  this function gives the non-relativistic KE on time")
        vSquared = np.multiply(self.velocity, self.velocity)
        vSquared = np.sum(vSquared, axis=-1)
        KE = 0.5*self.particle.m0*vSquared
        ax = plt.figure().add_subplot()
        ax.plot(self.time, vSquared)
        ax.set_title("Kinetic energy on time")

    def plotDeltaV(self):
        print("warning:  this function needs improvement")
        v = self.velocity
        deltav = np.zeros((np.shape(v)[0]-1))
        for i in range(1, len(v)-1):
            deltav[i] = np.linalg.norm(v[i] - v[i-1])
        
        ax = plt.figure().add_subplot()
        ax.plot(deltav)
        
        
        return
  
