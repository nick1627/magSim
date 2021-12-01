"""
This file contains particle objects, like protons and electrons
"""


class Particle:
    def __init__(self, mass, charge, position, velocity):
        self.m0 = mass          #Rest mass of particle
        self.q = charge         #charge of particle
        
        self.r = position
        self.v = velocity

    def getPosition(self):
        return self.r
    
    def getVelocity(self):
        return self.v

    def updatePosition(self, field, deltaT):
        #First need to get the next velocity
        #Then need to get the next position
        self.r = 1
        return

    

class Electron(Particle):
    def __init__(self, position, velocity):
        self.m0 = 9.11E-31
        self.q = -1.60E-19

        self.r = position
        self.v = velocity

class Proton(Particle):
    def __init__(self, position, velocity):
        self.m0 = 1.67E-27
        self.q = 1.60E-19

        self.r = position
        self.v = velocity