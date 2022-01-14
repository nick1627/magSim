"""
This file stores the code that analyses the simulation
"""
import numpy as np
from numpy.lib.arraysetops import isin
from fields import *
from particles import *
import time 

class Simulation:
    """
    This class manages everything to do with the simulation.  Simulations are
    comprised of a single particle moving within a magnetic field, and so a
    field and a particle object are taken as input.

    The timestep is the amount of time between successive particle position
    estimations.  Each simulation object has a particular timestep.

    Different runs of the simulation can have different durations, measured
    either in steps or in multiples of the characteristic timescale (the 
    gyro-period).
    """
    def __init__(self, field, particle, timestep, simDataPath = ""):
        #timestep always input as fraction of initial period of gyroradius

        #If previous simulation data not provided, we start a new simulation
        if simDataPath == "":
            self.timeStep = timestep
            self.field = field
            self.particle = particle

            #Diagnostic info
            self.position = []
            self.velocity = []
            self.time = []

            self.complete = False
        else:
            #We have prior data, so create simulation object based on that for analysis

            savedArrays = np.load(simDataPath)

            fieldArray = savedArrays["fieldData"]
            particleArray = savedArrays["particleData"]
            simulationArray = savedArrays["simulationData"]
            positionArray = savedArrays["positions"]
            velocityArray = savedArrays["velocities"]
            timeArray = savedArrays["times"]


            self.timeStep = simulationArray[0]
            self.field = "empty"
            self.particle = "empty"
            self.position = positionArray
            self.velocity = velocityArray
            self.time = timeArray
            self.complete = True
        
        return 
    
    def wipeSimData(self):
        self.position = []
        self.velocity = []
        self.time = []
        return


    def run(self, endTime=0, endStep=0, naturalUnits=True):
        """
        This function runs the simulation for a given duration in terms of time
        or steps.  Can run in either natural units or SI units.
        """
        #End time is always in units of time steps, not seconds
        
        if endTime == 0:
            if endStep == 0:
                raise(Exception("You didnt specify how the simulation should end!"))
            else:
                endOnTime = False
        else:
            endOnTime = True

        self.complete = False

        #first wipe the simulation data already existing.
        self.wipeSimData()

        print("Beginning simulation...")
        #Find start time
        startRealTime = time.time()
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

            if endOnTime:
                while currentTime < endTime: #
                    #increment time 
                    currentTime += self.timeStep
                    #update particle position
                    self.particle.updatePositionN(self.field, self.timeStep)
                    #Record data
                    self.position.append(self.particle.getPosition())
                    self.velocity.append(self.particle.getVelocity())
                    self.time.append(currentTime)


            else: #end based on step
                steps = 0
                while steps < endStep:
                    #increment time 
                    currentTime += self.timeStep
                    #update particle position
                    self.particle.updatePositionN(self.field, self.timeStep)
                    #Record data
                    self.position.append(self.particle.getPosition())
                    self.velocity.append(self.particle.getVelocity())
                    self.time.append(currentTime)

                    steps += 1
            
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

            if endOnTime == False:
                raise(Exception("Ending on step count has not been implemented yet for SI units"))

            #Record initial values 
            self.position.append(self.particle.getPosition())
            self.velocity.append(self.particle.getVelocity())
            self.time.append(currentTime)

            endTime *= self.particle.larmorPeriod
 
            if endOnTime:
                while currentTime < endTime: #
                    #increment time 
                    currentTime += self.timeStep*self.particle.larmorPeriod
                    #update particle position
                    self.particle.updatePositionSI(self.field, self.timeStep*self.particle.larmorPeriod)
                    #Record data
                    self.position.append(self.particle.getPosition())
                    self.velocity.append(self.particle.getVelocity())
                    self.time.append(currentTime)


            else: #end based on step
                steps = 0
                while steps < endStep:
                    #increment time 
                    currentTime += self.timeStep*self.particle.larmorPeriod
                    #update particle position
                    self.particle.updatePositionN(self.field, self.timeStep*self.particle.larmorPeriod)
                    #Record data
                    self.position.append(self.particle.getPosition())
                    self.velocity.append(self.particle.getVelocity())
                    self.time.append(currentTime)

                    steps += 1
            
            #Simulation has now completed
            self.position = np.array(self.position)
            self.velocity = np.array(self.velocity)

        elapsedRealTime = time.time() - startRealTime
        print("Simulation completed.  Time elapsed: " + str(elapsedRealTime) + " seconds.")

        self.complete = True

        return
        
    def saveData(self, filePath=""):
        """
        Saves all data required to recreate the simulation & plot graphs 
        """

        #First save data on field
        if isinstance(self.field, SHField):
            fieldArray = [self.field.a, self.field.g, self.field.h]
            fieldArray = np.array(fieldArray)
        elif isinstance(self.field, UniformField):
            fieldArray = np.array([self.field.B])
        else:
            fieldArray = np.array([])
        #Next save data on particle 
        initialPosition = self.position[0]
        initialVelocityDirection = self.velocity[0]/np.linalg.norm(self.velocity[0])

        gamma = 1/(np.sqrt(1 - (np.linalg.norm(self.velocity[0], axis=-1)/self.particle.c)**2))
        restMassEnergy = self.particle.m0*self.particle.c**2
        initialKE = (gamma - 1)*restMassEnergy #This is in joules

        particleArray = np.array([self.particle.m0, self.particle.q, initialPosition, initialVelocityDirection, initialKE])
        #Now save simulation output
        simulationArray = np.array([self.timeStep])

        # This saves the data 
        np.savez(filePath, fieldData = fieldArray, particleData = particleArray, simulationData = simulationArray, positions = self.position, velocities = self.velocity, times = self.time)
        return

    def plotLShellOnTime(self):
        #Need to get L from position
        L = np.zeros(np.shape(self.position)[0])
        for i in range(0, np.shape(self.position)[0]):
            sphericalPos = self.field.convertCartesianToPolar(self.position[i])
            L[i] = sphericalPos[0]/(np.sin(sphericalPos[1]))**2

        L = L/self.field.a
        #Now have obtained an array of L-shell
        #Now plot L-shell on time

        ax = plt.figure().add_subplot()
        ax.plot(self.time, L, color="purple")
        
        titleString = "L-shell on time"
        ax.set_title(titleString)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("L-shell (planetary radii)")

        return

    def plotFirstAIOnTime(self):
        #Plots the first adiabatic invariant on time
        r = self.position
        v = self.velocity
        m = self.particle.m0

        mu = np.zeros(np.shape(self.position)[0])

        for i in range(0, np.shape(r)[0]):
            B = self.field.getField(r[i])
            BMag = np.linalg.norm(B)
            BHat = B/BMag
            vPerp = v[i] - np.dot(v[i], BHat)*BHat

            mu[i] = (m/(2*BMag))*np.dot(vPerp, vPerp)
        
        ax = plt.figure().add_subplot()
        ax.plot(self.time, mu, color="blue")
        
        titleString = "1st Adiabatic Invariant on Time"
        ax.set_title(titleString)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("First Adiabatic Invariant")

        return

    def plotPositionOnTime(self, x=False, y=False, z=False):
        ax = plt.figure().add_subplot(projection = "3d")
        ax.plot(self.position[:,0], self.position[:,1], self.position[:,2])
        ax.set_title("Position of particle in 3D")

        if x:
            ax1 = plt.figure().add_subplot()
            ax1.plot(self.time, self.position[:, 0])
            ax1.set_title("x on time")
        if y:
            ax2 = plt.figure().add_subplot()
            ax2.plot(self.time, self.position[:, 1])
            ax2.set_title("y on time")
        if z:
            ax3 = plt.figure().add_subplot()
            ax3.plot(self.time, self.position[:, 2])
            ax3.set_title("z on time")

        ax4 = plt.figure().add_subplot()
        ax4.plot(self.position[:, 0], self.position[:, 1])
        ax4.set_aspect("equal")
        ax4.set_title("x and y positions as time evolves")

        return

    def plotKEOnTime(self):
        restMassEnergy = self.particle.m0*self.particle.c**2
        v = np.linalg.norm(self.velocity, axis=-1)
        gamma = 1/(np.sqrt(1 - (v/self.particle.c)**2))
        for g in gamma:
            if math.isnan(g) or math.isinf(g):
                print(g)
                raise(Exception("Warning!  Energy problems!"))

        Ek = (gamma - 1)*restMassEnergy #This is in joules
        Ek = Ek/(1000*1.6E-19)
        ax = plt.figure().add_subplot()
        ax.plot(self.time, Ek, color="red")
        gain = ((Ek[-1] - Ek[0])/Ek[0])*100
        gain = np.round(gain, decimals=3)
        titleString = "Kinetic Energy on Time - Percentage Gain: " + str(gain) + "%"
        ax.set_title(titleString)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Kinetic Energy (keV)")

        return

    def plotDeltaV(self):
        print("warning:  this function needs improvement")
        v = self.velocity
        deltav = np.zeros((np.shape(v)[0]-1))
        for i in range(1, len(v)-1):
            deltav[i] = np.linalg.norm(v[i] - v[i-1])
        
        ax = plt.figure().add_subplot()
        ax.plot(deltav)
        
        
        return

    def plotVelocityOnTime(self):
        ax = plt.figure().add_subplot()
        v = np.linalg.norm(self.velocity, axis =-1)
        ax.plot(self.time, v)
        gain = ((v[-1] - v[0])/v[0])*100
        gain = np.round(gain, decimals=3)
        titleString = "Speed on Time - Percentage Gain: " + str(gain) + "%"
        ax.set_title(titleString)

        return

    def plotVelocityErrorOnTime(self):
        ax = plt.figure().add_subplot()
        v = np.linalg.norm(self.velocity, axis=-1)
        v_error = abs(v - v[0])
        ax.plot(self.time, v_error)
        
        return
  
class SimulationManager:
    def __init__(self, fieldList, particleList, timeStepList, N = 10, mainFilePath = "Output/", fileNames = "auto", fileKeyWord = "", endTimeList = 0, endStepList = 0, naturalUnitsList = True):
        """
        fieldList:         A list of fields
        particleList:      A list of particles
        timeStepList:      A list of time steps to use, one for each sim
        N:                 The number of simulations.  Must match the lengths of the lists (or they may be extended)
        mainFilePath:      The directory you wish to save all the simulations to
        fileNames:         #The file names to be used.  Usually set to "auto" for automatic naming.  TODO: finish implementing  
        fileKeyWord:       A key word to appear in the file names to help you keep track of various different runs
        endTimeList:       A list of end times for the simulations, expressed in units of the timeStep
        endStepList:       A list of end step numbers for the sims
        naturalUnitsList:  A list of bool to determine whether natural units are to be used in the sims TODO: check
        """
        
        self.N = N

        if not isinstance(fieldList, list):
            fieldList = [fieldList]*self.N
        if not isinstance(particleList, list):
            particleList = [particleList]*self.N
        if not isinstance(timeStepList, list):
            timeStepList = [timeStepList]*self.N
        if not isinstance(endTimeList, list):
            endTimeList = [endTimeList]*self.N
        if not isinstance(endStepList, list):
            endStepList = [endStepList]*self.N
        if not isinstance(naturalUnitsList, list):
            naturalUnitsList = [naturalUnitsList]*self.N

        if len(fieldList) != self.N:
            raise(Exception("Length mismatch!"))
        if len(particleList) != self.N:
            raise(Exception("Length mismatch!"))
        if len(timeStepList) != self.N:
            raise(Exception("Length mismatch!"))
        if len(endTimeList) != self.N:
            raise(Exception("Length mismatch!"))
        if len(endStepList) != self.N:
            raise(Exception("Length mismatch!"))
        if len(naturalUnitsList) != self.N:
            raise(Exception("Length mismatch!"))


        self.fieldList = fieldList
        self.particleList = particleList
        self.timeStepList = timeStepList
        self.endTimeList = endTimeList
        self.endStepList = endStepList
        self.naturalUnitsList = naturalUnitsList

        #Now have lists of info set up for the simulations
        self.simulations = []
        for i in range(0, self.N):
            self.simulations.append(Simulation(fieldList[i], particleList[i], timeStepList[i]))

        #Deal with the filenames and paths
        self.filePaths = []
        if fileNames == "auto":
            for i in range(0, self.N):
                self.filePaths.append(mainFilePath + fileKeyWord + "-" + str(self.simulations[i].particle.name) + "-" + str(np.round(self.simulations[i].particle.initialEnergy)) + ".npz")

        return

    def runAllSims(self):
        print("Beginning simulations...")
        for i in range(0, self.N):
            print("Starting simulation %d of %d..." % (i+1, self.N))
            self.simulations[i].run(self.endTimeList[i], self.endStepList[i], self.naturalUnitsList[i])
            self.simulations[i].saveData(filePath = self.filePaths[i])

        print("All simulations complete.")

    def plotAllEnergy(self):
        print("Plotting particle energies...")
        for i in range(0, self.N):
            self.simulations[i].plotKEOnTime()





