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

    def plot3DField(self, xmin, xmax, ymin, ymax, zmin, zmax, noPoints, planetRadius, planetaryFilter = True, scaleOverride = 1):
        xstep = (xmax - xmin)/(noPoints - 1)
        ystep = (ymax - ymin)/(noPoints - 1)
        zstep = (zmax - zmin)/(noPoints - 1)

        x = np.zeros((noPoints, noPoints, noPoints))
        y = np.zeros(np.shape(x))
        z = np.zeros(np.shape(x))

        u = np.zeros(np.shape(x))
        v = np.zeros(np.shape(x))
        w = np.zeros(np.shape(x))


        biggestLength = 0
        initialPosition = np.array([xmin, ymin, zmin])
        #Define lattice vectors
        a = np.array([xstep, 0, 0])
        b = np.array([0, ystep, 0])
        c = np.array([0, 0, zstep])

        if planetaryFilter:
            for i in range(0, noPoints):
                for j in range(0, noPoints):
                    for k in range(0, noPoints):
                        r = initialPosition + i*a + j*b + k*c

                        x[i, j, k] = r[0]
                        y[i, j, k] = r[1]
                        z[i, j, k] = r[2]

                        if np.linalg.norm(r) > planetRadius:
                            B = self.getField(r)
                            currentLength = np.linalg.norm(B)
                            if currentLength > biggestLength:
                                biggestLength = currentLength
                        
                            u[i, j, k] = B[0]
                            v[i, j, k] = B[1]
                            w[i, j, k] = B[2]


        else:
            for i in range(0, noPoints):
                for j in range(0, noPoints):
                    for k in range(0, noPoints):
                        r = initialPosition + i*a + j*b + k*c
                      
                        x[i, j, k] = r[0]
                        y[i, j, k] = r[1]
                        z[i, j, k] = r[2]
                        
                        B = self.getField(r)
                        currentLength = np.linalg.norm(B)
                        if currentLength > biggestLength:
                            biggestLength = currentLength
                    
                        u[i, j, k] = B[0]
                        v[i, j, k] = B[1]
                        w[i, j, k] = B[2]


        scale = scaleOverride*xstep/biggestLength
        u = scale*u
        v = scale*v
        w = scale*w



        ax = plt.figure().add_subplot(projection = "3d")
        ax.quiver(x, y, z, u, v, w)
        #ax.set_aspect("equal")       

        return x, y, z, u, v, w
        

    def plot2DField(self, plane, intercept, noPoints, vecLength, planetRadius, vectors = "2D", scaleOverride = 1, planetaryFilter = True):
        #  This plots a 2D vector field (arrows may be in 3D) as a plane
        #  Normal vector of the plane is n, specified by user
        #  Plane will always intersect the origin

        #a is a vector in the plane
        #b is another vector in the plane, perpendicular to a
        #N is the number of points along one side of the plane

        if plane == "x":
            a = np.array([0, vecLength, 0])
            b = np.array([0, 0, vecLength])
            start = np.array([intercept, 0, 0])
        elif plane == "y":
            a = np.array([vecLength, 0, 0])
            b = np.array([0, 0, vecLength])
            start = np.array([0, intercept, 0])
        elif plane == "z":
            a = np.array([vecLength, 0, 0])
            b = np.array([0, vecLength, 0])
            start = np.array([0, 0, intercept])
        

        x = np.zeros((noPoints, noPoints))
        y = np.zeros(np.shape(x))
        z = np.zeros(np.shape(x))

        u = np.zeros(np.shape(x))
        v = np.zeros(np.shape(x))
        w = np.zeros(np.shape(x))

        start = start - ((noPoints-1)/2)*a - ((noPoints-1)/2)*b
        biggestLength = 0

        aHat = a/np.linalg.norm(a)
        bHat = b/np.linalg.norm(b)

        if planetaryFilter: #filter out vectors inside the planet
            if vectors == "2D":
                for i in range(0, np.shape(x)[0]):
                    for j in range(0, np.shape(x)[1]):
                        currentVector = start + i*a + j*b
                        if np.linalg.norm(currentVector) > planetRadius:

                            x[i, j] = currentVector[0]
                            y[i, j] = currentVector[1]
                            z[i, j] = currentVector[2]

                            currentB = self.getField(currentVector)
                            u[i, j] = currentB[0]
                            v[i, j] = currentB[1]
                            w[i, j] = currentB[2]

                            currentLength = pow(np.dot(currentB, aHat), 2) + pow(np.dot(currentB, bHat), 2)
                            if currentLength > biggestLength:
                                
                                biggestLength = currentLength
                
                biggestLength = np.sqrt(biggestLength)

            else: #vectors = 3D
                for i in range(0, np.shape(x)[0]):
                    for j in range(0, np.shape(x)[1]):
                        currentVector = start + i*a + j*b
                        if np.linalg.norm(currentVector) > planetRadius:

                            x[i, j] = currentVector[0]
                            y[i, j] = currentVector[1]
                            z[i, j] = currentVector[2]

                            currentB = self.getField(currentVector)
                            u[i, j] = currentB[0]
                            v[i, j] = currentB[1]
                            w[i, j] = currentB[2]

                            currentLength = np.linalg.norm(currentB)
                            if currentLength > biggestLength:
                                biggestLength = currentLength

        else: #Don't filter out vectors inside planet
            if vectors == "2D":
                for i in range(0, np.shape(x)[0]):
                    for j in range(0, np.shape(x)[1]):
                        currentVector = start + i*a + j*b
                        
                        x[i, j] = currentVector[0]
                        y[i, j] = currentVector[1]
                        z[i, j] = currentVector[2]

                        currentB = self.getField(currentVector)
                        u[i, j] = currentB[0]
                        v[i, j] = currentB[1]
                        w[i, j] = currentB[2]

                        currentLength = pow(np.dot(currentB, aHat), 2) + pow(np.dot(currentB, bHat), 2)
                        if currentLength > biggestLength:
                            biggestLength = currentLength
                
                biggestLength = np.sqrt(biggestLength)

            else:
                for i in range(0, np.shape(x)[0]):
                    for j in range(0, np.shape(x)[1]):
                        currentVector = start + i*a + j*b
                        
                        x[i, j] = currentVector[0]
                        y[i, j] = currentVector[1]
                        z[i, j] = currentVector[2]

                        currentB = self.getField(currentVector)
                        u[i, j] = currentB[0]
                        v[i, j] = currentB[1]
                        w[i, j] = currentB[2]
 
                        currentLength = np.linalg.norm(currentB)
                        if currentLength > biggestLength:
                            biggestLength = currentLength


        
        
        # u = scaleOverride*scale*u
        # v = scaleOverride*scale*v
        # w = scaleOverride*scale*w

        #Now we plot the vectors.  1 unit of B takes the lenght of biggestLength/vecLength metres (since x in metres)

        if vectors == "3D":
            ax2 = plt.figure().add_subplot(projection = "3d")
            ax2.quiver(x, y, z, u, v, w, scale = biggestLength/vecLength, scale_units = "x")
        else: #vectors = 2D
            ax2 = plt.figure().add_subplot()

            if plane == "x":
                ax2.quiver(y/planetRadius, z/planetRadius, v, w, scale = biggestLength*planetRadius/vecLength, scale_units = "x")
                
                
            elif plane == "y":
                ax2.quiver(x/planetRadius, z/planetRadius, u, w, scale = biggestLength*planetRadius/vecLength, scale_units = "x")
            elif plane == "z":
                ax2.quiver(x/planetRadius, y/planetRadius, u, v, scale = biggestLength*planetRadius/vecLength, scale_units = "x")

        ax2.set_aspect("equal")       

        return x, y, z, u, v, w



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

        self.nMax = np.shape(g)[0]
 

    def PnmCos(self, n, m, theta): #Returns Pnm(cos(theta)) for n up to 2
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
        
    def PnmCosDerivative(self, n, m, theta): #Returns derivative of Pnm(cos(theta)) wrt theta for n up to 2
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
              
                Br += (n+1)*frac2*(self.g[n-1, m]*np.cos(m*phi) + self.h[n-1, m]*np.sin(m*phi))*self.PnmCos(n, m, theta)
                Btheta -= (frac2*(self.g[n-1, m]*np.cos(m*phi) + self.h[n-1, m]*np.sin(m*phi))*self.PnmCosDerivative(n, m, theta))
                Bphi += m*frac2*(self.g[n-1, m]*np.sin(m*phi) - self.h[n-1, m]*np.cos(m*phi))*self.PnmCos(n, m, theta)
        
        Bphi *= 1/np.sin(theta)

        #Now need to recover cartesian mag field components
        Bspherical = np.array([Br, Btheta, Bphi])
        T = np.array([[np.sin(theta)*np.cos(phi),   np.cos(theta)*np.cos(phi),  -np.sin(theta)], 
                        [np.sin(theta)*np.sin(phi), np.cos(theta)*np.sin(phi),  np.cos(theta)],
                        [np.cos(theta),             -np.sin(theta),             0]])

        Bcartesian = np.matmul(T, Bspherical)

        return np.array(Bcartesian)


class UniformField(Field):
    def __init__(self, fieldVec):
        self.B = fieldVec
    
    def getField(self):
        return self.B

