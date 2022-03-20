"""
This file stores the code that analyses the simulation
"""

import sys, os
sys.path.insert(0, os.getcwd())

import numpy as np
from fields import *
from particles import *
import time 
import scipy as sp
import tools

#TODO:  make initial energy a property of the simulation,not the particle

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
    def __init__(self, field=0, particle=0, stepsPerPeriod=50, simDataPath = "", initialPhase=-1):
        #timestep always input as fraction of initial period of gyroradius

        #If previous simulation data not provided, we start a new simulation
        if simDataPath == "":
            
            self.stepsPerPeriod: int = stepsPerPeriod
            self.field = field
            self.particle = particle
            self.initialPhase = initialPhase #stored but never used.  Negative values imply phase was never given

            #Diagnostic info
            #Position, velocity and time may be entered into the array in natural units, but will
            #always eventually get converted to SI units.
            self.position = []
            self.velocity = []
            self.time = []

            self.complete = False
        else:
            
            #We have prior data, so create simulation object based on that for analysis

            savedArrays = np.load(simDataPath, allow_pickle=True)

            fieldArray = savedArrays["fieldData"]
            particleArray = savedArrays["particleData"]
            simulationArray = savedArrays["simulationData"]
            positionArray = savedArrays["positions"]
            velocityArray = savedArrays["velocities"]
            timeArray = savedArrays["times"]


            self.timeStep = simulationArray[0]
            self.initialPhase = simulationArray[1]
            if np.shape(fieldArray)[0] == 4:
                #then we should have spherical harmonic field
                self.field = SHField(fieldArray[0], fieldArray[1], fieldArray[2], 0, 0, fieldArray[3])
            else:
                raise Exception("This has not been dealt with.  Something aobut uniform fields")
                self.field = "empty"
            
            if particleArray[0] == sp.constants.m_p:
                #proton
                self.particle = Proton(particleArray[2], particleArray[3], particleArray[4])
            elif particleArray[0] == sp.constants.m_e:
                #electron
                self.particle = Electron(particleArray[2], particleArray[3], particleArray[4])
            else:
                self.particle = Particle(particleArray[0], particleArray[1], particleArray[2], particleArray[3], particleArray[4])
            
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

    def getLarmorPeriod(self):
        #get magnitude of B at position of particle
        Bmag = np.linalg.norm(self.field.getField(self.particle.getPosition()))
        #get velocity vector of particle
        v = self.particle.getVelocity() #returns velocity in SI units
        gamma = 1/np.sqrt(1 - (np.linalg.norm(v)/sp.constants.c)**2)
        larmorPeriod = gamma*self.particle.m0*2*np.pi/(abs(self.particle.q)*Bmag)
        return larmorPeriod

    def getLarmorRadius(self, index, alternativeEnergy=""):
        """
        Computes the larmor radius at a point in the simulation after the simulation has run.

        index:  integer
        alternativeEnergy:  eV
        """
        if self.complete == False:
            raise(Exception("Simulation must be complete for this function to run"))

        if self.field.rotationFlag == "R":
            self.field.rotate("Field")

        if alternativeEnergy=="":
            position = self.position[index]
            velocity = self.velocity[index]

         
            B = self.field.getField(position)
            Bmag = np.linalg.norm(B)
            gamma = 1/np.sqrt(1-(np.linalg.norm(velocity)/sp.constants.c)**2)
            Bdir = B/Bmag
            v_perp = velocity - np.dot(velocity, Bdir)*Bdir
            larmorRadius = (gamma*self.particle.m0*np.linalg.norm(v_perp))/(abs(self.particle.q)*Bmag)
        
        else:
            #we're given an alternative energy
            position=self.position[index]
            velocityDirection = self.velocity[index]/np.linalg.norm(self.velocity[index])

            m0 = self.particle.m0
            
            newE = alternativeEnergy*sp.constants.e #newE is now in joules
            newVelocity = (np.sqrt(1-((self.particle.m0*sp.constants.c**2)/(self.particle.m0*sp.constants.c**2 + newE))**2))*sp.constants.c*velocityDirection

            B = self.field.getField(position)
            Bmag = np.linalg.norm(B)
            gamma = 1/np.sqrt(1-(np.linalg.norm(newVelocity)/sp.constants.c)**2)
            Bdir = B/Bmag
            v_perp = newVelocity - np.dot(newVelocity, Bdir)*Bdir
            larmorRadius = (gamma*self.particle.m0*np.linalg.norm(v_perp))/(abs(self.particle.q)*Bmag)

        return larmorRadius


    def getChangeInLarmorRadius(self, index, alternativeEnergy):
        """
        require alternative energy in eV
        """
        actualLarmorRadius = self.getLarmorRadius(index)
        newLarmorRadius = self.getLarmorRadius(index, alternativeEnergy=alternativeEnergy)

        return newLarmorRadius - actualLarmorRadius
    
    def run(self, endStep=1000):
        """
        This is a more presumptuous way of running the simulation.  It does not include
        options for having SI units, as natural units are used as standard.  Furthermore,
        there is no option for setting the in-simulation end time; instead only steps are
        used to measure when the simulation should finish.

        The timestep is adaptive.  Every stepsPerPeriod, the method updates the current
        gyro-period and sets the timestep such that there are still stepsPerPeriod steps
        in said period.
        """

        self.complete = False
        self.wipeSimData()

        print("Beginning simulation...")
        initialRunTime = time.time()

        currentTime = 0
        currentStep = 0

        #Store the initial values of time, position, speed
        self.position.append(self.particle.getPosition())
        self.velocity.append(self.particle.getVelocity(True))
        self.time.append(currentTime)

        #currentTimeStep is the time in seconds each step is worth.
        currentLarmorPeriod = self.getLarmorPeriod()
        currentTimeStep = currentLarmorPeriod/self.stepsPerPeriod
        while currentStep < endStep:
            currentStep += 1
            #Branch depending on whether the current step is a multiple of 50 (or stepsPerPeriod)
            if currentStep % self.stepsPerPeriod == 0:
                #Need to update the time step
                currentLarmorPeriod = self.getLarmorPeriod()
                currentTimeStep = currentLarmorPeriod/self.stepsPerPeriod
            
            #increment time
            currentTime += currentTimeStep #currentTimeStep should update automatically because of python being weird
   
            Bdash = self.field.getField(self.particle.r)*(self.particle.q*currentLarmorPeriod/self.particle.m0)
            #update particle position using differential equation solutions
            self.particle.updatePositionN(Bdash, 1/self.stepsPerPeriod, currentLarmorPeriod)
            #Record data
            self.position.append(self.particle.getPosition())
            self.velocity.append(self.particle.getVelocity(True))
            self.time.append(currentTime)



        elapsedRunTime = time.time() - initialRunTime
        print("Simulation completed.  Time elapsed: " + str(elapsedRunTime) + " seconds.")

        self.complete = True


        #Finalise the data
        #Convert to numpy arrays for easy maths
        self.position = np.array(self.position)
        self.velocity = np.array(self.velocity)
        self.time = np.array(self.time)
        #Now ensure stored in SI
        self.velocity = self.velocity*sp.constants.c

        return

        
    def saveData(self, filePath=""):
        """
        Saves all data required to recreate the simulation & plot graphs 
        """

        #First save data on field
        if isinstance(self.field, SHField):
            fieldArray = [self.field.a, self.field.g, self.field.h, self.field.getDipoleFlag()]
            fieldArray = np.array(fieldArray, dtype=object)
        elif isinstance(self.field, UniformField):
            fieldArray = np.array([self.field.B])
        else:
            fieldArray = np.array([])
        #Next save data on particle 
        initialPosition = self.position[0]
        initialVelocityDirection = self.velocity[0]/np.linalg.norm(self.velocity[0])

        gamma = 1/(np.sqrt(1 - (np.linalg.norm(self.velocity[0], axis=-1)/sp.constants.c)**2))
        restMassEnergy = self.particle.m0*sp.constants.c**2
        initialKE = (gamma - 1)*restMassEnergy #This is in joules
        initialKE = initialKE/sp.constants.e

        particleArray = np.array([self.particle.m0, self.particle.q, initialPosition, initialVelocityDirection, initialKE], dtype=object)
        #Now save simulation output
        simulationArray = np.array([self.stepsPerPeriod, self.initialPhase])

        # This saves the data 
        np.savez(filePath, fieldData = fieldArray, particleData = particleArray, simulationData = simulationArray, positions = self.position, velocities = self.velocity, times = self.time)
        return

    def printSpecs(self, detailLevel = 1):
        print("The specifications for this simulation are:")
        
        return

    def plotLShellOnTime(self, titleAddition="", otherSims = [], legendList=[], finalIndex = -1):

        if otherSims == []:
            #Need to get L from position
            L = np.zeros(np.shape(self.position[:finalIndex])[0])
            for i in range(0, np.shape(L)[0]):
                sphericalPos = self.field.convertCartesianToPolar(self.position[i])
                L[i] = sphericalPos[0]/(np.sin(sphericalPos[1]))**2

            L = L/self.field.a
            #Now have obtained an array of L-shell
            #Now plot L-shell on time

            ax = plt.figure().add_subplot()
            ax.plot(self.time[:finalIndex], L, color="purple")
            
            titleString = "L-shell on time" + " " + titleAddition
            ax.set_title(titleString)
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("L-shell (planetary radii)")

        else: #we have multiple simulations
            sims = [self] + otherSims
            L = np.zeros((len(sims), np.shape(self.position)[0]))

            minTime = sims[0].time[-1]
            for i in range(1, len(sims)):
                if sims[i].time[-1] < minTime:
                    minTime = sims[i].time[-1]

            for i in range(0, len(sims)):
                for j in range(0, np.shape(self.position)[0]):
                    sphericalPos = self.field.convertCartesianToPolar(sims[i].position[j])
                    L[i, j] = sphericalPos[0]/(np.sin(sphericalPos[1]))**2

            L = L/self.field.a

            ax = plt.figure().add_subplot()
            colourList = ["blue", "purple", "red", "green"]
            for i in range(0, len(sims)):
                ax.plot(sims[i].time, L[i, :], color=colourList[i], label = legendList[i])
            
            titleString = "L-shell on time" + " " + titleAddition

            plt.xlim((0, minTime))
            ax.set_title(titleString)
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("L-shell (planetary radii)")
            ax.legend()

        return

    def plotFirstAIOnTime(self, titleAddition = ""):
        #Plots the first adiabatic invariant on time
        r = self.position
        v = self.velocity
        m = self.particle.m0

        if isinstance(self.field, SHField):
            self.field.rotate("Field")

        mu = np.zeros(np.shape(self.position)[0])

        for i in range(0, np.shape(r)[0]):
            B = self.field.getField(r[i])
            # print(self.field.nMax)
            BMag = np.linalg.norm(B)
            BHat = B/BMag
            vPerp = v[i] - np.dot(v[i], BHat)*BHat
            gammaSquared = 1/(1-(np.linalg.norm(v[i])/sp.constants.c)**2)

            mu[i] = (m/(2*BMag))*np.dot(vPerp, vPerp)*gammaSquared
        
        ax = plt.figure().add_subplot()
        ax.plot(self.time, mu, color="blue")
        
        titleString = "1st Adiabatic Invariant on Time" + " " + titleAddition
        ax.set_title(titleString)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("First Adiabatic Invariant")

        return

    def plotPositionOnTime(self, x=False, y=False, z=False):
        ax = plt.figure().add_subplot(projection = "3d")
        ax.plot(self.position[:,0], self.position[:,1], self.position[:,2], color="red")

        #REMOVE THESE LINES IF U WANT
        scale = np.linalg.norm(self.position[0,:])
        ax.plot([0,0], [0,0], [-scale, scale], color="w", linestyle="none")
        ax.plot([0,0], [-scale, scale], [0,0], color="w", linestyle="none")
        ax.plot([-scale, scale], [0,0], [0,0], color="w", linestyle="none")


        if isinstance(self.field, SHField):
            a = self.field.a
            # draw sphere
            u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
            x1 = a*np.cos(u)*np.sin(v)
            y1 = a*np.sin(u)*np.sin(v)
            z1 = a*np.cos(v)
            ax.plot_wireframe(x1, y1, z1, color="blue")
        ax.set_title("Position of particle over time")

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

    def plotKEOnTime(self, titleAddition="", plot=True, returnData=False):
        restMassEnergy = self.particle.m0*sp.constants.c**2
        v = np.linalg.norm(self.velocity, axis=-1)
        gamma = 1/(np.sqrt(1 - (v/sp.constants.c)**2))
        for g in gamma:
            if math.isnan(g) or math.isinf(g):
                print(g)
                raise(Exception("Warning!  Energy problems!"))

        Ek = (gamma - 1)*restMassEnergy #This is in joules
        Ek = Ek/(1000*sp.constants.e)
        ax = plt.figure().add_subplot()

        if plot:
            ax.plot(self.time, Ek, color="red")
            gain = ((Ek[-1] - Ek[0])/Ek[0])*100
            gain = np.round(gain, decimals=3)
            titleString = "Kinetic Energy on Time - Percentage Gain: " + str(gain) + "%" + " " + titleAddition
            ax.set_title(titleString, y=1.04)
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("Kinetic Energy (keV)")

        if returnData:
            return Ek*1000 #return KE in eV
        else:
            return

    def getKE(self, index, unit="eV"):
        """
        Returns the relativistic KE at a given index.
        """
        if self.complete == False:
            raise(Exception("Simulation must be complete for this function to work"))
        restMassEnergy = self.particle.m0*sp.constants.c**2
        v = np.linalg.norm(self.velocity[index], axis=-1)
        gamma = 1/np.sqrt(1 - (v/sp.constants.c)**2)
        Ek = (gamma-1)*restMassEnergy

        if unit=="SI":
            return Ek
        elif unit=="eV":
            return Ek/sp.constants.e
        else:
            raise(Exception("Unrecognised unit!"))


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

    def getFirstEquatorialIndex(self):
        """
        Returns the first index after the particle has crossed the magnetic equator
        """
        if self.complete==False:
            raise(Exception("Simulation has not yet finished.  Simulation must be complete before this function can run."))

        if abs(self.position[0, 2]) > 10**(-5):
            # print(self.position[0])
            raise(Exception("Particle does not start at magnetic equator, so cannot perform analysis."))

        #first collect up all the data
        #Index of initial event is 0
        #Need to find final index

        found = False
        counter=0
        length = np.shape(self.position)[0]       
        if self.position[1, 2] > 0:
            # print("Particle heading north")
            #particle heading North initially
            while (found == False) and (counter < length):
                if self.position[counter, 2]<0:
                    #particle has crossed back into the other hemisphere
                    found = True
                else:
                    counter+=1
        else:
            # print("Particle heading south")
            #particle heading South initially
            counter = 5
            while (found == False) and (counter < length):
                if self.position[counter, 2]>0:
                    #particle has crossed back into other hemisphere
                    found = True
                else:
                    counter+=1
        
        if found == False:
            raise(Exception("The particle in this simulation did not cross back into the other hemisphere.  No results can be found."))

        # print("The index is " + str(counter))
        if counter < 100:
            print("yeah it's here")
            print(counter)
            print("this is a bit small right?")

        #counter is now the index of the final point.

        return counter


    def saveBounceData(self, filePath):
        """
        This function can run after the simulation has completed.  It saves data from the first bounce in a given file.
        """
        if self.complete==False:
            raise(Exception("Simulation has not yet finished.  Simulation must be complete before this function can run."))

        if abs(self.position[0, 2]) > 10**(-5):
            # print(self.position[0])
            raise(Exception("Particle does not start at magnetic equator, so cannot perform analysis."))

        #first collect up all the data
        #Index of initial event is 0
        #Need to find final index
        counter = self.getFirstEquatorialIndex()


        #counter is now the index of the final point.
                
            
        initialKE = self.getKE(0, unit="eV")
        initialGyroradius = self.getLarmorRadius(0)

        finalGyroradius = self.getLarmorRadius(counter)

        initialRadius = np.linalg.norm(self.getGuidingCentrePosition(0, initialGyroradius))
        finalPos = self.getGuidingCentrePosition(counter, finalGyroradius)
        # print(finalPos)
        finalRadius = np.linalg.norm(finalPos)

        positionError, velocityError = self.getPositionAndVelocityError(counter)

        finalKE = self.getKE(counter, unit="eV")
     
        #calculate pitch angle
        v = self.velocity[0]
        d = np.sqrt(v[0]**2 + v[1]**2)
        pitchAngle = np.arctan((d/v[2]))

        initial_gc = self.getGuidingCentrePosition(0, initialGyroradius)
        initialPhi = np.arctan2(initial_gc[1], initial_gc[0])

        tools.saveRegionData2(filePath, "N", self.particle.name, not(self.field.dipoleOnly), initialKE, finalKE, pitchAngle, self.initialPhase, initialRadius, finalRadius, initialGyroradius, finalGyroradius, positionError, initialPhi)

        return

    def getGuidingCentrePosition(self, index, larmorRadius):
        """
        Computes the position of the guiding centre
        """
        if self.complete == False:
            raise(Exception("The simulation must be complete"))

        perpDirection = self.particle.q*np.cross(self.velocity[index], self.field.getField(self.position[index]))
        perpDirection = perpDirection/np.linalg.norm(perpDirection)
        guidingCentrePosition = self.position[index] + larmorRadius*perpDirection
        return guidingCentrePosition

    def plotAltitudeOnTime(self):

        
        alts = np.linalg.norm(self.position, axis=1)

        
        ax1 = plt.figure().add_subplot()
        ax1.plot(self.time, alts)
        ax1.set_title("altitude (m) on time")
       
        return       

    # def getEquatorialGuidingCentreError(self):
    #     """
    #     This function finds the error in the position of the guiding centre when the particle first crosses the equator
    #     """

    #     equatorialIndex = self.getFirstEquatorialIndex()
    #     firstIndex = equatorialIndex - 50
    #     secondIndex = equatorialIndex + 50

    #     try:
    #         temp = self.position[secondIndex]
    #     except:
    #         raise(Exception("Not enough points to do the analysis.  Modify it"))

    #     closePositions = self.position[firstIndex:secondIndex, :]
    #     guidingCentres = np.zeros(np.shape(closePositions))
    #     for i in range(0, np.shape(guidingCentres)[0]):
    #         rL = self.getLarmorRadius(firstIndex + i)
    #         guidingCentres[i, :] = self.getGuidingCentrePosition(firstIndex + i, larmorRadius=rL)


    #     # plt.figure()
    #     # plt.plot(closePositions[:, 0], closePositions[:, 1], color="blue")
    #     # plt.plot(guidingCentres[:, 0], guidingCentres[:, 1], color="red")
    #     # plt.show()    
    #     ax = plt.figure().add_subplot(projection = "3d")
    #     ax.plot(closePositions[:,0], closePositions[:,1], closePositions[:,2])
    #     ax.plot(guidingCentres[:, 0], guidingCentres[:, 1], guidingCentres[:, 2])
    #     plt.show()
    #     #Now we need to find the variation in phi for these positions in the spherical coord system
    #     for i in range(0, np.shape(guidingCentres)[0]):
    #         guidingCentres[i,:] = self.convertCartesianToPolar(guidingCentres[i, :])
        
    #     #now the positions are of the form r, theta, phi
    #     #check phi values go from 0 to 2pi
    #     guidingCentres[:,2] = guidingCentres[:,2] % 2*np.pi
    #     #We're looking for the max and min phi values 
    #     maxPhi = np.max(guidingCentres[:, 2])
    #     minPhi = np.min(guidingCentres[:, 2])
    #     deltaPhi = abs(maxPhi - minPhi)
    #     if deltaPhi > np.pi:
    #         deltaPhi = 2*np.pi - deltaPhi
    #     deltaPhi = abs(deltaPhi) #just in case
    #     # print(deltaPhi)
    #     #now we turn this delta phi into a length
    #     rL = self.getLarmorRadius(equatorialIndex)
    #     gcPosition = self.getGuidingCentrePosition(equatorialIndex, rL)
    #     gcRadius = np.linalg.norm(gcPosition)
    #     print(gcRadius/25600000)
    #     gcError = gcRadius*np.tan(deltaPhi/2)
    #     gcError = 2*gcError

    #     return gcError

    def getPositionAndVelocityError(self, N):
        """
        Returns the position and velocity error at the given index.

        Must be done after the simulation has run.

        N:      The index we are calculating the error for
        """
        #first check the simulation has run
        if self.complete == False:
            raise(Exception("The simulation must be complete"))

     

        #find timestep to use
        timeSteps = np.zeros(N-1)
        for i in range(0, N-1):
            timeSteps[i] = self.time[i+1] - self.time[i]

        h = np.max(timeSteps)    

        #Need to find when the maximum acceleration was happening, which hopefully will be at the pole
       
            
        # averageAcceleration = accelerationSum/N
        # accelerations = np.zeros(N)
        # for i in range(0, N):
        #     accelerations[i] = np.linalg.norm((self.velocity[i+1] - self.velocity[i])/(self.time[i+1] - self.time[i]))
           
        # averageAcceleration = np.mean(accelerations)

        #now maxTimeStep is the maximum time step that occurred.
        #Define various computational errors
        mu_x0 = 10**(-7)
        mu_v0 = 10**(-18)
        mu_x = 10**(-7)
        mu_v = 10**(-18)

        #Define the truncation errors
        a_v = h**7
        # a_x = 0.5*averageAcceleration*h**2

        velocityError = (mu_v0 + (N-1)*(a_v + mu_v))*sp.constants.c


        # positionError = mu_x0 + (N-1)*(a_x + mu_x) + h*sp.constants.c*(mu_x0 + a_v*(N-1)*(N-2)/2 - a_v*N + mu_v*(N-1)*(N-2)/2 - mu_v*N)
        accelerationTerm = 0
        for i in range(0, N):
            current_h = self.time[i+1] - self.time[i]
            current_a = (self.velocity[i+1] - self.velocity[i])/(self.time[i+1] - self.time[i])
            accelerationTerm += 0.5*current_a*current_h**2

        accelerationTerm = np.array(accelerationTerm[0], accelerationTerm[1])
        accelerationTerm = np.linalg.norm(accelerationTerm)

        # print("terms")
        # print(mu_x0)
        # print((N-1)*mu_x)
        # print(accelerationTerm)
        # print(h*sp.constants.c*mu_x0)
        # print(h*sp.constants.c*(a_v*(N-1)*(N-2)/2 - a_v*N))
        # print(h*sp.constants.c*(mu_v*(N-1)*(N-2)/2 - mu_v*N))
        # print("terms end")

        # print("Comparing dumb truncation errors: x, then v:")
        # print((N-1)*a_x)
        # print(h*sp.constants.c*(a_v*(N-1)*(N-2)/2 - a_v*N))
        # print("comparison over")

        # positionError = mu_x0 + accelerationTerm + (N-1)*(mu_x) + h*sp.constants.c*(mu_x0 + a_v*(N-1)*(N-2)/2 - a_v*N + mu_v*(N-1)*(N-2)/2 - mu_v*N)
        # positionError = mu_x0 + accelerationTerm + (N-1)*(mu_x) + h*sp.constants.c*(mu_x0 + mu_v*(N-1)*(N-2)/2 - mu_v*N)
        positionError = accelerationTerm

        return positionError, velocityError


    def convertCartesianToPolar(self, rvec, theta=0, phi=0, origin=True):
        if origin:
        
            r = np.sqrt(np.power(rvec[0], 2) + np.power(rvec[1], 2) + np.power(rvec[2], 2))
            theta = np.arctan2(np.sqrt(np.power(rvec[0], 2) + np.power(rvec[1], 2)), rvec[2])
            phi = np.arctan2(rvec[1], rvec[0])

            result = np.array([r, theta, phi])
        
        else:
            Sp = np.sin(phi)
            St = np.sin(theta)
            Cp = np.cos(phi)
            Ct = np.cos(theta)

            T = np.array([[St*Cp, St*Sp, Ct], 
                            [Ct*Cp, Ct*Sp, -St],
                            [-Sp, Cp, 0]])

            result = np.matmul(T, rvec)
        
        return result


    def convertPolarToCartesian(self, rvec, theta, phi):
        #rvec is the components of the vector to be converted in order (a, b, c) for a*rhat + b*theta hat + c*phihat
        #For vectors from the origin, rvec should read (r, 0, 0), but theta and phi must be given too

        T = np.array([[np.sin(theta)*np.cos(phi),   np.cos(theta)*np.cos(phi),  -np.sin(phi)], 
                        [np.sin(theta)*np.sin(phi), np.cos(theta)*np.sin(phi),  np.cos(phi)],
                        [np.cos(theta),             -np.sin(theta),             0]])

        return np.matmul(T, rvec)

    def convertCartesianToPolar(self, rvec, theta=0, phi=0, origin=True):
        if origin:
        
            r = np.sqrt(np.power(rvec[0], 2) + np.power(rvec[1], 2) + np.power(rvec[2], 2))
            theta = np.arctan2(np.sqrt(np.power(rvec[0], 2) + np.power(rvec[1], 2)), rvec[2])
            phi = np.arctan2(rvec[1], rvec[0])

            result = np.array([r, theta, phi])
        
        else:
            Sp = np.sin(phi)
            St = np.sin(theta)
            Cp = np.cos(phi)
            Ct = np.cos(theta)

            T = np.array([[St*Cp, St*Sp, Ct], 
                            [Ct*Cp, Ct*Sp, -St],
                            [-Sp, Cp, 0]])

            result = np.matmul(T, rvec)
        
        return result


    def getFirstBounceThetaPhi(self):

        maxIndex = self.getFirstEquatorialIndex()

        positions = self.position[:maxIndex, :]
        # print(positions[4,:])

        phiTheta = np.zeros((np.shape(positions)[0], 2))
        for i in range(0, np.shape(positions)[0]):
            # print(positions[i, :])
            newvec = self.convertCartesianToPolar(positions[i,:])
            phiTheta[i, :] = np.array([newvec[2], newvec[1]]) * (180/np.pi)
            phiTheta[i, 0] = phiTheta[i, 0] % 360
            # phiTheta[i, 1] = phiTheta[i, 1] % 180
            # print(phiTheta[i, :])

        # print(phiTheta[:, 1])
        
        return phiTheta
  
class SimulationManager:
    """
    A simulation manager mangages multiple simulations
    """
    def __init__(self, fieldList, particleList, stepsPerPeriodList=50, N = 10, mainFilePath = "Output/", fileKeyWord = "", endStepList = 0, initialPhaseList=-1):
        """
        fieldList:          A list of fields
        particleList:       A list of particles
        stepsPerPeriodList: A list of the number of steps per period to use, one for each sim.
        N:                  The number of simulations.  Must match the lengths of the lists (or they may be extended)
        mainFilePath:       The directory you wish to save all the simulations to
        fileNames:          #The file names to be used.  Usually set to "auto" for automatic naming.  TODO: finish implementing  
        fileKeyWord:        A key word to appear in the file names to help you keep track of various different runs
        endStepList:        A list of end step numbers for the sims
        """
        
        self.N = N

        if not isinstance(fieldList, list):
            fieldList = [fieldList]*self.N
        if not isinstance(particleList, list):
            particleList = [particleList]*self.N
        if not isinstance(stepsPerPeriodList, list):
            stepsPerPeriodList = [stepsPerPeriodList]*self.N
        if not isinstance(endStepList, list):
            endStepList = [endStepList]*self.N
        if not isinstance(initialPhaseList, list):
            initialPhaseList = [initialPhaseList]*self.N
   

        if len(fieldList) != self.N:
            raise(Exception("Length mismatch!"))
        if len(particleList) != self.N:
            raise(Exception("Length mismatch!"))
        if len(stepsPerPeriodList) != self.N:
            raise(Exception("Length mismatch!"))
        if len(endStepList) != self.N:
            raise(Exception("Length mismatch!"))
        if len(initialPhaseList) != self.N:
            raise(Exception("Length mismatch!"))



        self.fieldList = fieldList
        self.particleList = particleList
        self.stepsPerPeriodList = stepsPerPeriodList
        self.endStepList = endStepList
        self.initialPhaseList = initialPhaseList

        #Now have lists of info set up for the simulations
        self.simulations = []
        for i in range(0, self.N):
            self.simulations.append(Simulation(fieldList[i], particleList[i], stepsPerPeriodList[i], initialPhase=initialPhaseList[i]))

        #Deal with the filenames and paths
        self.filePaths = []
     
        for i in range(0, self.N):
            particleName = str(self.simulations[i].particle.name)
            initEnergy = str(np.round(self.simulations[i].particle.initialEnergy))
            if isinstance(self.fieldList[i], UniformField):
                fieldType="uniformField"
            elif isinstance(self.fieldList[i], SHField):
                if self.fieldList[i].dipoleOnly:
                    fieldType="dipoleOnly"
                else:
                    fieldType="fullField"
            else:
                fieldType=""

            self.filePaths.append(mainFilePath + fileKeyWord + "-" + particleName + "-" + fieldType + "-" + initEnergy + ".npz")

        return

        

    def runAllSims(self):
        print("Beginning simulations...")
        for i in range(0, self.N):
            print("Starting simulation %d of %d..." % (i+1, self.N))
            self.simulations[i].run(self.endStepList[i])
            self.simulations[i].saveData(filePath = self.filePaths[i])

        print("All simulations complete.")
        return

    def plotAllEnergy(self):
        print("Plotting particle energies...")
        for i in range(0, self.N):
            self.simulations[i].plotKEOnTime()
        return

class LocationCheck(SimulationManager):
    def __init__(self, L, theta, phi, N, energyList, particleType, gyroPhase, field, endStepList=500000, fileNameAddition=""):
        """
        This performs a set of simulations where the particle is aimed at a particular location

        L:              The initial L shell of the particles
        theta:          Theta of target location (degrees)
        phi:            Phi of target location (degrees)
        N:              Number of particles sent in/simulation runs
        energyList:     List of energies to use for each run
        particleType:   The type of particle ("proton" or "electron")
        gyroPhase:      The phase of the particle in radians
        field:          The field object.  The field must be an SHField.
        """

        fileNameAddition += (str(np.round(L)) + "-" + str(np.round(theta)) + "-" + str(np.round(phi)) + "-" + str(np.round(gyroPhase)) + "-")

        particleList = []


        centre_r = L*field.a
        centre_theta = (np.pi/180)*90
        centre_phi = (np.pi/180)*copy.copy(phi)
        centre_x = centre_r*np.sin(centre_theta)*np.cos(centre_phi)
        centre_y = centre_r*np.sin(centre_theta)*np.sin(centre_phi)
        centre_z = centre_r*np.cos(centre_theta)
        guidingCentrePosition = np.array([centre_x, centre_y, centre_z])

        phi = (np.pi/180)*phi

        initialB = np.linalg.norm(field.getField(guidingCentrePosition))
        
        fieldStatus = copy.copy(field.getDipoleFlag())
        field.setDipoleOnly(True)
        tempB = field.getField(guidingCentrePosition)
        
        field.setDipoleOnly(fieldStatus)

        latitude = (np.pi/180)*(90 - theta)
        gyroPhase = (np.pi/180)*gyroPhase
        
        gyroPhase = gyroPhase % (2*np.pi)
        

        # #calculate the equatorial pitch angle alpha
        alpha = np.arcsin(np.sqrt((np.cos(latitude)**6)/np.sqrt(1 + 3*(np.sin(latitude))**2)))


        if particleType == "proton":
            #need to create the protons
            for i in range(0, N):
                #Need to do different calcs for each, since they have different energies
                m0 = sp.constants.m_p
                q = sp.constants.e

                Ek = energyList[i]*sp.constants.e #Ek is now in joules
                speedSI = (np.sqrt(1-((m0*sp.constants.c**2)/(m0*sp.constants.c**2 + Ek))**2))*sp.constants.c


                # #compute larmor radius
                larmorRadius = self.getLarmorRadius(initialB, speedSI, alpha, m0, q)

                #Now set the direction modifier depending on the field direction and the particle charge
                if tempB[2]>0:
                    directionModifier = -np.pi/2 #because of the sign of the particle charge
                else:
                    directionModifier = np.pi/2
                
                velocityDirection = np.array([np.sin(alpha)*np.cos(phi + gyroPhase + directionModifier), np.sin(alpha)*np.sin(phi + gyroPhase + directionModifier), np.sign(latitude)*np.cos(alpha)])
                position = guidingCentrePosition + larmorRadius*np.array([np.cos(phi + gyroPhase), np.sin(phi + gyroPhase), 0])   

                particleList.append(Proton(position, velocityDirection, energyList[i])) 

        elif particleType == "electron":
            #Need to create list of electrons
            for i in range(0, N):

                #Need to do different calcs for each, since they have different energies
                m0 = sp.constants.m_e
                q = -sp.constants.e

                Ek = energyList[i]*sp.constants.e #Ek is now in joules
                speedSI = (np.sqrt(1-((m0*sp.constants.c**2)/(m0*sp.constants.c**2 + Ek))**2))*sp.constants.c


                # #compute larmor radius
                larmorRadius = self.getLarmorRadius(initialB, speedSI, alpha, m0, q)

                #Now set the direction modifier depending on the field direction and the particle charge
                if tempB[2] > 0:
                    directionModifier = np.pi/2 #because of the sign of the particle charge
                else:
                    directionModifier = -np.pi/2
                
                velocityDirection = np.array([np.sin(alpha)*np.cos(phi + gyroPhase + directionModifier), np.sin(alpha)*np.sin(phi + gyroPhase + directionModifier), np.sign(latitude)*np.cos(alpha)])
                position = guidingCentrePosition + larmorRadius*np.array([np.cos(phi + gyroPhase), np.sin(phi + gyroPhase), 0])   

                particleList.append(Electron(position, velocityDirection, energyList[i]))
        else:
            raise(Exception("Invalid particle type"))

       
        #400 steps definitely works for e- at 1MeV
        super(LocationCheck, self).__init__(field, particleList, stepsPerPeriodList=50, N=N, fileKeyWord="locationCheck" + fileNameAddition, endStepList=endStepList, initialPhaseList=gyroPhase)


    def getLarmorRadius(self, B, v, alpha, m0, q):
        #B in teslas, vector
        #v scalar in m/s
        #alpha is pitch angle, in radians

        
        
        v = abs(v)
       
        v_perp = v*np.sin(alpha)
        gamma = 1/np.sqrt(1 - (v/sp.constants.c)**2)
        larmorRadius = (gamma*m0*v_perp)/(abs(q)*B)
        larmorRadius = abs(larmorRadius)
        return larmorRadius
