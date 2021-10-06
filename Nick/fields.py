"""
Module for storing magnetic field code.

There is an overarching class called Field, from which SHField is derived (along with other fields)


"""

import numpy as np
import pyshtools as sh

class Field:
    # def __init__(self, cartesian = True):
    #     #basic field superclass to deal with coordinate system conversions
    #     #we assume that the field wants to be in cartesian at the base level
    
    def polarToCartesianB(self, B):
        #input B in the form Br, Btheta, Bphi in a vector B
        #output Bx, By, Bz
        Bx = B[0]*np.sin(B[1])*np.cos(B[2])
        By = B[0]*np.sin(B[1])*np.sin(B[2])
        Bz = B[0]*np.cos(B[1])

        return np.array([Bx, By, Bz])

    def cartesianToPolarB(self, B):

        Br = np.sqrt(B[0]**2 + B[1]**2 + B[2]**2)
        Btheta = np.arctan2(np.sqrt(B[0]**2 + B[1]**2), B[2])
        Bphi = np.arctan2(B[1], B[0])

        return np.array([Br, Btheta, Bphi])




class SHField(Field):
    def __init__(self, radius, g, h, G, H, g_error, h_error, G_error, H_error):
        #This is a planetary magnetic field using the spherical harmonic method
        self.a = radius

        self.g = g
        self.h = h
        self.G = G
        self.H = H

        self.nMax_i = 1 #Need to get this from the shape of g etc
        self.nMax_e = 1

    
    def getField(self, rvec, coordSys = "spherical"):
        #Assume spherical coordinate system
        #position is input as r, theta, phi

        r = rvec[0]
        theta = rvec[1]
        phi = rvec[2]

        #Analytical form given by Connerney (1993)

        #Will return B as a vector containing Br, Btheta, Bphi

        Br = 0
        Btheta = 0
        Bphi = 0

        frac = self.a/r

        for n in range(1, self.nMax_i+1):
            frac2 = frac**(n+2)
            for m in range(0, n+1):
                
                Br += (n+1)*frac2*(self.g[n, m]*np.cos(m*phi) + self.h[n,m]*np.sin(m*phi))*sh.PlmSchmidt(np.cos(theta))
                Btheta -= (frac2*(self.g[n, m]*np.cos(m*phi) + self.h[n,m]*np.sin(m*phi))*DERIVATIVE)
                Bphi += m*frac2*(self.g[n, m]*np.sin(m*phi) - self.h[n,m]*np.cos(m*phi))*sh.PlmSchmidt(np.cos(theta))
        
        Bphi *= 1/np.sin(theta)


        return np.array([Br, Btheta, Bphi])


class UniformField(Field):
    def __init__(self, fieldVec):
        self.B = fieldVec
    
    def getField(self):
        return self.B

