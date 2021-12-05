"""
This file contains particle objects, like protons and electrons
"""

from os import times
from typing import AsyncContextManager
import numpy as np
import copy
import math
class Particle:
    def __init__(self, mass, charge, position, velocity):
        self.m0 = mass          #Rest mass of particle
        self.q = charge         #charge of particle
        
        self.r = position
        self.v = velocity       #velocity of particle in natural units (divide by c)

        self.naturalUnits = False
        self.c = 3E8
        return

    def computeConversionFactors(self, B):
        #No natural units used yet
        #B in TESLAS
        Bmag = np.linalg.norm(B)
        Bdir = B/Bmag
        v_perp = self.v - np.dot(self.v, Bdir)*Bdir
        # if self.naturalUnits:
        #     v_perp = c*v_perp

        self.larmorRadius = (self.m0*np.linalg.norm(v_perp))/(abs(self.q)*Bmag)
        self.larmorPeriod = self.m0/(abs(self.q)*Bmag)
        return

    def convertToNatural(self):
        #Assume everything in SI, and need to convert to natural units
        self.r = self.r/self.larmorRadius
        self.v = self.v/self.c
        self.naturalUnits = True
       

    def convertToSI(self):
        self.r = self.r*self.larmorRadius
        self.v = self.v*self.c
        self.naturalUnits = False

    def getPosition(self):
        return self.r
    
    def getVelocity(self):
        return self.v

    def updatePositionN(self, field, timeStep):
        #First need to get the next velocity
        vdash = copy.copy(self.v)
        Bdash = field.getField(self.r)
        #Find k values
        root21 = np.sqrt(21)
        k = np.zeros((8,3))
        k[1, :] = self.accelerationN(vdash, Bdash) #k1
        k[2, :] = self.accelerationN(vdash + timeStep*k[1, :], Bdash) #k2
        k[3, :] = self.accelerationN(vdash + timeStep*(0.125*(3*k[1, :] + k[2, :])), Bdash) #k3
        k[4, :] = self.accelerationN(vdash + timeStep*((8*k[1, :] + 2*k[2, :] + 8*k[3, :])/27), Bdash) #k4
        k[5, :] = self.accelerationN(vdash + timeStep*(((9*root21 - 21)*k[1, :] - 8*(7 - root21)*k[2, :] + 48*(7 - root21)*k[3, :] - 3*(21 - root21)*k[4, :])/392), Bdash) #k5
        k[6, :] = self.accelerationN(vdash + timeStep*((-5*(231 + 51*root21)*k[1, :] - 40*(7 + root21)*k[2, :] - 320*root21*k[3, :] + 3*(21 + 121*root21)*k[4, :] + 392*(6+root21)*k[5, :])/1960), Bdash) #k6
        k[7, :] = self.accelerationN(vdash + timeStep*((15*(22+7*root21)*k[1, :] + 120*k[2, :] + 40*(7*root21 - 5)*k[3, :] - 63*(3*root21 - 2)*k[4, :] - 14*(49 + 9*root21)*k[5, :] + 70*(7-root21)*k[6, :])/180), Bdash)
        #dont need k7
        
        self.v = vdash + 0.2*timeStep*((16/27)*k[1, :] + (6656/2565)*k[3, :] + (28561/11286)*k[4, :] - 0.9*k[5, :] + (2/11)*k[6, :])#(9*k[1] + 64*k[3] + 49*k[5] + 49*k[6])
        #self.v = self.v + timeStep*(9*k[1] + 64*k[3] + 49*k[5] + 49*k[6] + 9*k[7])/180
        # self.v = vdash + timeStep*(self.acceleration(vdash, Bdash))
        #Then need to get the next position
   

        # print(np.linalg.norm(self.v))
        
        self.r = self.r + (self.larmorPeriod*3E8)*self.v*timeStep/self.larmorRadius  #0.2*timeStep*((16/27)*k[1] + (6656/2565)*k[3] + (28561/11286)*k[4] - 0.9*k[5] + (2/11)*k[6])
        print(np.linalg.norm(self.r))
        return

    def updatePositionSI(self, field, timeStep):
        #First need to get the next velocity
        v = copy.copy(self.v)
        B = field.getField(self.r)
        #Find k values
        root21 = np.sqrt(21)
        k = np.zeros((8,3))
        k[1, :] = self.accelerationSI(v, B) #k1
        k[2, :] = self.accelerationSI(v + timeStep*k[1, :], B) #k2
        k[3, :] = self.accelerationSI(v + timeStep*(0.125*(3*k[1, :] + k[2, :])), B) #k3
        k[4, :] = self.accelerationSI(v + timeStep*((8*k[1, :] + 2*k[2, :] + 8*k[3, :])/27), B) #k4
        k[5, :] = self.accelerationSI(v + timeStep*(((9*root21 - 21)*k[1, :] - 8*(7 - root21)*k[2, :] + 48*(7 - root21)*k[3, :] - 3*(21 - root21)*k[4, :])/392), B) #k5
        k[6, :] = self.accelerationSI(v + timeStep*((-5*(231 + 51*root21)*k[1, :] - 40*(7 + root21)*k[2, :] - 320*root21*k[3, :] + 3*(21 + 121*root21)*k[4, :] + 392*(6+root21)*k[5, :])/1960), B) #k6
        k[7, :] = self.accelerationSI(v + timeStep*((15*(22+7*root21)*k[1, :] + 120*k[2, :] + 40*(7*root21 - 5)*k[3, :] - 63*(3*root21 - 2)*k[4, :] - 14*(49 + 9*root21)*k[5, :] + 70*(7-root21)*k[6, :])/180), B)
        #dont need k7
        
        self.v = v + 0.2*timeStep*((16/27)*k[1, :] + (6656/2565)*k[3, :] + (28561/11286)*k[4, :] - 0.9*k[5, :] + (2/11)*k[6, :])
        # self.v = v + timeStep*self.accelerationSI()
       
        #Then need to get the next position
        self.r = self.r + self.v*timeStep

        return

        
    def accelerationN(self, vdash, Bdash):
        return np.cross(vdash, Bdash)*np.sqrt(1 - np.linalg.norm(vdash)**2)
    
    def accelerationSI(self, v, B):
        return (self.q/self.m0)*np.sqrt(1 - np.linalg.norm(v/self.c)**2)*np.cross(v, B)

    # def speed(self, vdash):
    #     return self.larmorPeriod*3E8*vdash/self.larmorRadius
    

class Electron(Particle):
    def __init__(self, position, velocity):
        self.m0 = 9.11E-31
        self.q = -1.60E-19

        self.r = position
        self.v = velocity

        self.c = 3E8

        return

class Proton(Particle):
    def __init__(self, position, velocity):
        self.m0 = 1.67E-27
        self.q = 1.60E-19

        self.r = position
        self.v = velocity

        self.c = 3E8

        return