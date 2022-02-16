"""
This file contains particle objects, like protons and electrons
"""

from os import times
import numpy as np
import copy
import math
import scipy as sp

#TODO:  make initial energy a property of the simulation,not the particle
class Particle:
    """
    This class is a general class for a particle.  To create a new particle, input
    the mass, charge, position, velocity direction vector and kinetic energy of the 
    particle.
    """
    def __init__(self, mass, charge, position, velocityDirection=np.array([0, 0, 1]), kineticEnergy=1000000, particleName = "None"):
        """
        mass:               Mass of particle in kg
        charge:             Charge of particle in C
        position:           Position of particle as a vector in m
        velocityDirection:  The direction that the particle should travel in initially, expressed as a cartesian vector (usually)
        kineticEnergy:      The initial KE of the particle in eV
        particleName:       e.g. "proton" or "neutron"
        targetSetup:        Boolean.  True means you're trying to aim for a particular position in the magnetosphere, in which case the
                            position and velocityDirection variables will mean something else.
        """
        #kinetic energy in eV
        self.m0 = mass          #Rest mass of particle in kg
        self.q = charge         #charge of particle in coulombs

        self.r = position       #position in metres

        #normalise velocity direction just in case
        velocityDirection = velocityDirection/np.linalg.norm(velocityDirection)
        

        #velocity is stored as v/c ONLY
        Ek = kineticEnergy*sp.constants.e #Ek is now in joules
        self.v = (np.sqrt(1-((self.m0*sp.constants.c**2)/(self.m0*sp.constants.c**2 + Ek))**2))*velocityDirection

        self.name = particleName
        self.initialEnergy = kineticEnergy #Stored in eV


        return


    def getLarmorRadius(self, B, v, alpha):
        #B in teslas, vector
        #v scalar in m/s
        #alpha is pitch angle, in radians

        print(B, v, alpha)
        
        v = abs(v)
       
        v_perp = v*np.sin(alpha)
        gamma = 1/np.sqrt(1 - (v/sp.constants.c)**2)
        larmorRadius = (gamma*self.m0*v_perp)/(abs(self.q)*B)
        larmorRadius = abs(larmorRadius)
        return larmorRadius


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
        self.r = copy.copy(self.r) + (larmorPeriod*sp.constants.c)*self.v*timeStep#*np.sqrt(1 - np.linalg.norm(self.v)**2)
        
        return


       
    def accelerationN(self, vdash, Bdash):
        return np.cross(vdash, Bdash)*np.sqrt(1 - np.linalg.norm(vdash)**2)



class Electron(Particle):
    """
    A subclass of particle, with the mass and charge already set to 
    that of an electron.
    """
    def __init__(self, position, velocityDirection, kineticEnergy):

        super(Electron, self).__init__(sp.constants.m_e, -sp.constants.e, position, velocityDirection, kineticEnergy, particleName="Electron")




        return

class Proton(Particle):
    """
    A subclass of particle, with the mass and charge already set to
    that of a proton.
    """
    def __init__(self, position, velocityDirection, kineticEnergy):

        super(Proton, self).__init__(sp.constants.m_p, sp.constants.e, position, velocityDirection, kineticEnergy, particleName="Proton")


        return