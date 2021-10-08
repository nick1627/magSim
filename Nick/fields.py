"""
Module for storing magnetic field code.

There is an overarching class called Field, from which SHField is derived (along with other fields)


"""

import numpy as np
import pyshtools as sh
import matplotlib.pyplot as plt

class Field:
    def __init__(self):
        self.B = 0
        
        #     #basic field superclass to deal with coordinate system conversions
        #     #we assume that the field wants to be in cartesian at the base level
        return
    
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

    def getField(self, posVec):
        #This is a placeholder
        return self.B

    def plotField(self, xmin, xmax, ymin, ymax, zmin, zmax, xstep = 1000, ystep = 1000, zstep = 1000):
        x = np.arange(xmin, xmax + xstep, xstep)
        y = np.arange(ymin, ymax + ystep, ystep)
        z = np.arange(zmin, zmax + zstep, zstep)

        grid_x, grid_y, grid_z = np.meshgrid(x, y, z)
        u, v, w = np.array(np.shape(grid_x))

        for i in range(0, np.shape(grid_x)[0]):
            for j in range(0, np.shape(grid_x)[1]):
                for k in range(0, np.shape(grid_x)[2]):
                    r = np.array([grid_x[i, j, k], grid_y[i, j, k], grid_z[i, j, k]])
                    B = self.getField(r)
                    u[i, j, k] = B[0]
                    v[i, j, k] = B[1]
                    w[i, j, k] = B[2]



        # for xval in x:
        #     for yval in y:
        #         for zval in z:
        #             r = np.array([xval, yval, zval])
        #             [u, v, w] = self.getField(r) 


        ax = plt.figure().add_subplot(projection = "3d")
        ax.quiver(x, y, z, u, v, w)
        




class SHField(Field):
    def __init__(self, radius, g, h, g_error, h_error):
        #This is a planetary magnetic field using the spherical harmonic method
        self.a = radius

        self.g = g
        self.h = h
        # self.G = G
        # self.H = H

        self.g_error = g_error
        self.h_error = h_error

        self.nMax = 1 #Need to get this from the shape of g etc
 

    def PnmCos(n, m, theta): #Returns Pnm(cos(theta)) for n up to 2
        if n == 2:
            if m == 0:
                return 1.5*(pow(np.cos(theta), 2) - 1/3)
            elif m == 1:
                return np.sqrt(3)*np.cos(theta)*np.sin(theta)
            elif m == 2:
                return 0.5*np.sqrt(3)*pow(np.sin(theta), 2)
            else:
                raise Exception("Error:  m = " + m + " is invalid!")

        elif n == 1:
            if m == 0:
                return np.cos(theta)
            elif m == 1:
                return np.sin(theta)
            else:
                raise Exception("Error:  m = " + m + " is invalid!")

        else:
            raise Exception("Error:  n = " + n + " is invalid!")
        
    def PnmCosDerivative(n, m, theta): #Returns derivative of Pnm(cos(theta)) wrt theta for n up to 2
        if n == 2:
            if m == 0:
                return -3*np.sin(theta)*np.cos(theta)
            elif m == 1:
                return np.sqrt(3)*np.cos(2*theta)
            elif m == 2:
                return np.sqrt(3)*np.sin(theta)*np.cos(theta)
            else:
                raise Exception("Error:  m = " + m + " is invalid!")

        elif n == 1:
            if m == 0:
                return -np.sin(theta)
            elif m == 1:
                return np.cos(theta)
            else:
                raise Exception("Error:  m = " + m + " is invalid!")

        else:
            raise Exception("Error:  n = " + n + " is invalid!")


    
    def getField(self, rvec):
        #Input rvec will always be cartesian


        #Assume spherical coordinate system
        #we want positions in r, theta, phi 

        r = np.sqrt(pow(rvec[0], 2) + pow(rvec[1], 2) + pow(rvec[2], 2))
        theta = np.arctan2(np.sqrt(pow(rvec[0], 2) + pow(rvec[1], 2)), rvec[2])
        phi = np.arctan2(rvec[1], rvec[0])

        #Analytical form given by Connerney (1993)

        #Will get B as a vector containing Br, Btheta, Bphi

        Br = 0
        Btheta = 0
        Bphi = 0

        frac = self.a/r

        for n in range(1, self.nMax+1):
            frac2 = frac**(n+2)
            for m in range(0, n+1):
                #g[n,m] difference in matrix index to n value
                print(n, m, theta)
                Br += (n+1)*frac2*(self.g[n-1, m]*np.cos(m*phi) + self.h[n-1,m]*np.sin(m*phi))*self.PnmCos(n, m, theta)
                Btheta -= (frac2*(self.g[n-1, m]*np.cos(m*phi) + self.h[n-1,m]*np.sin(m*phi))*self.PnmCosDerivative(n, m, theta))
                Bphi += m*frac2*(self.g[n-1, m]*np.sin(m*phi) - self.h[n-1,m]*np.cos(m*phi))*self.PnmCos(n, m, theta)
        
        Bphi *= 1/np.sin(theta)

        #Now need to recover cartesian mag field components
        Bspherical = np.array([Br, Btheta, Bphi])
        T = np.array([[np.sin(theta)*np.cos(phi), np.cos(theta)*np.cos(phi), -np.sin(theta)], 
                        [np.sin(theta)*np.sin(phi), np.cos(theta)*np.sin(phi), np.cos(theta)],
                        [np.cos(theta), -np.sin(theta), 0]])

        Bcartesian = np.matmul(T, Bspherical)

        return np.array(Bcartesian)


class UniformField(Field):
    def __init__(self, fieldVec):
        self.B = fieldVec
    
    def getField(self):
        return self.B

