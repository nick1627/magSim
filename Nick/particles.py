"""
This file contains particle objects, like protons and electrons
"""

from os import times
import numpy as np
import copy
import math
import scipy as sp
class Particle:
    """
    This class is a general class for a particle.  To create a new particle, input
    the mass, charge, position, velocity direction vector and kinetic energy of the 
    particle.
    """
    def __init__(self, mass, charge, position, velocityDirection, kineticEnergy, particleName = "None"):
        #kinetic energy in eV
        self.m0 = mass          #Rest mass of particle in kg
        self.q = charge         #charge of particle in coulombs
        
        self.r = position       #position in metres

        # self.naturalUnits = False

        velocityDirection = velocityDirection/np.linalg.norm(velocityDirection)
        Ek = kineticEnergy*sp.constants.e #Ek is now in joules

        #velocity is stored as v/c ONLY
        self.v = (np.sqrt(1-((self.m0*sp.constants.c**2)/(self.m0*sp.constants.c**2 + Ek))**2))*velocityDirection

        self.name = particleName
        self.initialEnergy = kineticEnergy #Stored in eV


        return

    # def computeConversionFactors(self, B):
    #     #No natural units used yet
    #     #B in TESLAS
    #     Bmag = np.linalg.norm(B)
    #     Bdir = B/Bmag
    #     v_perp = self.v - np.dot(self.v, Bdir)*Bdir
    #     gamma = 1/np.sqrt(1 - (np.linalg.norm(self.v))**2)
       
    #     self.larmorRadius = (gamma*self.m0*np.linalg.norm(v_perp))/(abs(self.q)*Bmag)
    #     self.larmorPeriod = gamma*self.m0*2*np.pi/(abs(self.q)*Bmag)

    #     return

    # def convertToNatural(self):
    #     #Assume everything in SI, and need to convert to natural units
    #     self.r = self.r/self.larmorRadius
    #     self.v = self.v/sp.constants.c
    #     self.naturalUnits = True
       

    # def convertToSI(self):
    #     self.r = self.r*self.larmorRadius
    #     self.v = self.v*sp.constants.c
    #     self.naturalUnits = False

    def getPosition(self):
        return self.r
    
    def getVelocity(self, natural=False):
        if natural:
            return self.v
        else:
            return sp.constants.c*self.v

    def updatePositionN(self, Bdash, timeStep, larmorPeriod):
        #Bdash is the dimensionless value of the B field at the current location of the particle
        #timeStep is the step size into the future that is being made, as a fraction of the larmor period.
        #larmorPeriod in seconds
        
        #First need to get the next velocity
        #Find k values
        root21 = np.sqrt(21)
        k = np.zeros((8,3))
        k[1, :] = self.accelerationN(self.v, Bdash) #k1
        k[2, :] = self.accelerationN(self.v + timeStep*k[1, :], Bdash) #k2
        k[3, :] = self.accelerationN(self.v + timeStep*(0.125*(3*k[1, :] + k[2, :])), Bdash) #k3
        k[4, :] = self.accelerationN(self.v + timeStep*((8*k[1, :] + 2*k[2, :] + 8*k[3, :])/27), Bdash) #k4
        k[5, :] = self.accelerationN(self.v + timeStep*(((9*root21 - 21)*k[1, :] - 8*(7 - root21)*k[2, :] + 48*(7 - root21)*k[3, :] - 3*(21 - root21)*k[4, :])/392), Bdash) #k5
        k[6, :] = self.accelerationN(self.v + timeStep*((-5*(231 + 51*root21)*k[1, :] - 40*(7 + root21)*k[2, :] - 320*root21*k[3, :] + 3*(21 + 121*root21)*k[4, :] + 392*(6+root21)*k[5, :])/1960), Bdash) #k6
        k[7, :] = self.accelerationN(self.v + timeStep*((15*(22+7*root21)*k[1, :] + 120*k[2, :] + 40*(7*root21 - 5)*k[3, :] - 63*(3*root21 - 2)*k[4, :] - 14*(49 + 9*root21)*k[5, :] + 70*(7-root21)*k[6, :])/180), Bdash)
        
        #update velocity.  Note that this v is actually v/c
        self.v = self.v + timeStep*(9*k[1, :] + 64*k[3, :] + 49*k[5, :] + 49*k[6, :] + 9*k[7, :])/180
 
        #Then update position.  This r is r in metres.  
        self.r = copy.copy(self.r) + (larmorPeriod*sp.constants.c)*self.v*timeStep*np.sqrt(1 - np.linalg.norm(self.v)**2)
        
        return


       
    def accelerationN(self, vdash, Bdash):
        return np.cross(vdash, Bdash)*np.sqrt(1 - np.linalg.norm(vdash)**2)


    # def updatePositionSI(self, field, timeStep):
    #     #First need to get the next velocity
    #     v = copy.copy(self.v)
    #     B = field.getField(self.r)
    #     #Find k values
    #     root21 = np.sqrt(21)
    #     k = np.zeros((8,3))
    #     k[1, :] = self.accelerationSI(v, B) #k1
    #     k[2, :] = self.accelerationSI(v + timeStep*k[1, :], B) #k2
    #     k[3, :] = self.accelerationSI(v + timeStep*(0.125*(3*k[1, :] + k[2, :])), B) #k3
    #     k[4, :] = self.accelerationSI(v + timeStep*((8*k[1, :] + 2*k[2, :] + 8*k[3, :])/27), B) #k4
    #     k[5, :] = self.accelerationSI(v + timeStep*(((9*root21 - 21)*k[1, :] - 8*(7 - root21)*k[2, :] + 48*(7 - root21)*k[3, :] - 3*(21 - root21)*k[4, :])/392), B) #k5
    #     k[6, :] = self.accelerationSI(v + timeStep*((-5*(231 + 51*root21)*k[1, :] - 40*(7 + root21)*k[2, :] - 320*root21*k[3, :] + 3*(21 + 121*root21)*k[4, :] + 392*(6+root21)*k[5, :])/1960), B) #k6
    #     k[7, :] = self.accelerationSI(v + timeStep*((15*(22+7*root21)*k[1, :] + 120*k[2, :] + 40*(7*root21 - 5)*k[3, :] - 63*(3*root21 - 2)*k[4, :] - 14*(49 + 9*root21)*k[5, :] + 70*(7-root21)*k[6, :])/180), B)
    #     #dont need k7
        
    #     self.v = v + 0.2*timeStep*((16/27)*k[1, :] + (6656/2565)*k[3, :] + (28561/11286)*k[4, :] - 0.9*k[5, :] + (2/11)*k[6, :])
    #     # self.v = v + timeStep*self.accelerationSI()
       
    #     #Then need to get the next position
    #     self.r = self.r + self.v*timeStep

    #     return

 
    
    # def accelerationSI(self, v, B):
    #     return (self.q/self.m0)*np.sqrt(1 - np.linalg.norm(v/sp.constants.c)**2)*np.cross(v, B)


    

class Electron(Particle):
    """
    A subclass of particle, with the mass and charge already set to 
    that of an electron.
    """
    def __init__(self, position, velocityDirection, kineticEnergy):
        #kinetic energy in keV
        self.m0 = sp.constants.m_e
        self.q = -sp.constants.e

        self.r = position
    

        velocityDirection = velocityDirection/np.linalg.norm(velocityDirection)
        Ek = kineticEnergy*sp.constants.e #Ek is now in joules
        self.v = (np.sqrt(1-((self.m0*sp.constants.c**2)/(self.m0*sp.constants.c**2 + Ek))**2))*velocityDirection

        self.name = "Electron"
        self.initialEnergy = kineticEnergy #in eV

        return

class Proton(Particle):
    """
    A subclass of particle, with the mass and charge already set to
    that of a proton.
    """
    def __init__(self, position, velocityDirection, kineticEnergy):
        self.m0 = sp.constants.m_p
        self.q = sp.constants.e

        self.r = position #in metres
    

        velocityDirection = velocityDirection/np.linalg.norm(velocityDirection)
        Ek = kineticEnergy*sp.constants.e #Ek is now in joules
        self.v = (np.sqrt(1-((self.m0*sp.constants.c**2)/(self.m0*sp.constants.c**2 + Ek))**2))*velocityDirection

        self.name = "Proton"
        self.initialEnergy = kineticEnergy #in eV

        return