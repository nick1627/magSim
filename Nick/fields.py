"""
Module for storing magnetic field code.

There is an overarching class called Field, from which SHField is derived (along with other fields)

TODO:  make sure code can handle field at pole somehow

TODO:  Convert all graph-plotting code to the following format:
        1.  All start with plot rather than get data
        2.  An array is input to define which bits to plot, eg an array of longitudes or L-shells
        3.  For each element of the array, the data will be found and the graph is optionally plotted
        4.  The data is returned

"""

import numpy as np
# from numpy.lib.arraysetops import isin
#import pyshtools as sh
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LogNorm
import copy

class Field:
    def __init__(self):
        self.B = 0
        self.naturalUnits = False
        
        #     #basic field superclass to deal with coordinate system conversions
        #     #we assume that the field wants to be in cartesian at the base level
        return


    

    def getUnitState(self):
        return self.naturalUnits

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


    def rotate2D(self, vec, theta):
        R = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
        return(np.matmul(R, vec))


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

        ax = plt.figure().add_subplot(projection = "3d")
        ax.quiver(x, y, z, u*scale, v*scale, w*scale)
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

        # Now we plot the vectors.  1 unit of B takes the lenght of biggestLength/vecLength metres (since x in metres)

        if vectors == "3D":
            ax2 = plt.figure().add_subplot(projection = "3d")
            ax2.quiver(x, y, z, u, v, w, scale = biggestLength/vecLength, scale_units = "x")
        else: #vectors = 2D
            ax2 = plt.figure().add_subplot()

            if plane == "x":
                ax2.quiver(y/planetRadius, z/planetRadius, v, w, scale = biggestLength*planetRadius/vecLength, scale_units = "x")
                ax2.set_xlabel("y/a")
                ax2.set_ylabel("z/a")
                ax2.set_title("Plane perpendicular to x-axis")
            elif plane == "y":
                ax2.quiver(x/planetRadius, z/planetRadius, u, w, scale = biggestLength*planetRadius/vecLength, scale_units = "x")
                ax2.set_xlabel("x/a")
                ax2.set_ylabel("z/a")
                ax2.set_title("Plane perpendicular to y-axis")
            elif plane == "z":
                ax2.quiver(x/planetRadius, y/planetRadius, u, v, scale = biggestLength*planetRadius/vecLength, scale_units = "x")
                ax2.set_xlabel("x/a")
                ax2.set_ylabel("y/a")
                ax2.set_title("Plane perpendicular to z-axis")

        ax2.set_aspect("equal")       

        return x, y, z, u, v, w



    



class SHField(Field):
    def __init__(self, radius, g, h, g_error=0, h_error=0, dipoleOnly=False):
        #This is a planetary magnetic field using the spherical harmonic method
        self.a = radius

        self.g = g
        self.h = h
        # self.G = G
        # self.H = H

        self.g_error = g_error
        self.h_error = h_error

        self.dipoleOnly = dipoleOnly
        if self.dipoleOnly:
            self.nMax = 1
        else:
            self.nMax = np.shape(g)[0]

        self.true_nMax = np.shape(g)[0]

        self.rotationFlag = "R"
        self.R = np.identity(3)
        self.Rinv = np.identity(3)

        self.naturalUnits = False

    def setDipoleOnly(self, dipoleFlag):
        """
        Set the dipoleFlag to true to have a field with only dipole components
        """
        if dipoleFlag == True:
            self.nMax = 1
        else:
            self.nMax = copy.deepcopy(self.true_nMax)

        self.dipoleOnly = dipoleFlag

        return

    def getDipoleFlag(self):
        return self.dipoleOnly

    def setNaturalUnits(self, activate, charge, restMass, period): #I believe this should be removed...
        #This switches the units to being natural units
        if activate:
            self.g = (charge*period/restMass)*self.g
            self.h = (charge*period/restMass)*self.h
            self.naturalUnits = True
        else:
            self.g = (restMass/(charge*period))*self.g
            self.h = (restMass/(charge*period))*self.h
            self.naturalUnits = False
        return    

    def rotate(self, rotationKey): #Sets the rotation matrix for the whole system.
        #rotationKey is a string that tells you what system to rotate into.
        #rotationKey can take values "Field" and "Rotation", or "F" and "R"

        if rotationKey == ("Field" or "F"):
            self.rotationFlag = "F"
        elif rotationKey == ("Rotation" or "R"):
            self.rotationFlag = "R"
        else:
            raise(Exception("Warning!  Invalid rotation set."))
        

        if self.rotationFlag == "R":
            self.R = np.identity(3)
            self.Rinv = np.identity(3)
        else: #In this case self.rotationFlag = "F", field-aligned coordinates.

            #This uses the components of g and h to align the dipole of the field with the z axis
            #First get the vector of the dipole moment
            #g10 = z, g11 = x, h11 = y
            zF = np.array([self.g[0, 1], self.h[0, 1], self.g[0, 0]]) #z axis in new system
            zF = zF/np.linalg.norm(zF) #should now be unit vector
            
            zR = np.array([0, 0, 1])

            xF = np.cross(zF, zR)
            xF = xF/np.linalg.norm(xF)

            yF = np.cross(zF, xF)
            yF = yF/np.linalg.norm(yF)


            self.R = np.transpose(np.vstack((xF, yF, zF)))


            self.Rinv = np.linalg.inv(self.R)

        return
 

  

    def PnmCos(self, n, m, theta): #Returns Pnm(cos(theta)) for n up to 2
        if n == 2:
            if m == 0:
                return 1.5*(np.power(np.cos(theta), 2) - 1/3)
            elif m == 1:
                return np.sqrt(3)*np.cos(theta)*np.sin(theta)
            elif m == 2:
                return 0.5*np.sqrt(3)*np.power(np.sin(theta), 2)
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
                return 0.5*np.sqrt(3)*np.sin(2*theta)
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


    def PnmCosDerivative2(self, n, m, theta): #Returns second derivative of Pnm(cos(theta)) wrt theta for n up to 2
        if n == 2:
            if m == 0:
                return -3*np.cos(2*theta)
            elif m == 1:
                return -2*np.sqrt(3)*np.sin(2*theta)
            elif m == 2:
                return np.sqrt(3)*np.cos(2*theta)
            else:
                raise Exception("Error:  m = " + m + " is invalid!")

        elif n == 1:
            if m == 0:
                return -np.cos(theta)
            elif m == 1:
                return -np.sin(theta)
            else:
                raise Exception("Error:  m = " + m + " is invalid!")

        else:
            raise Exception("Error:  n = " + n + " is invalid!")
        

    def getField(self, rvec, returnCartesian=True, returnR = False): #Returns the cartesian magnetic field at the specified cartesian position (either frame)
        #Input rvec will always be cartesian
        #Input can either be in F frame or R frame.  
        #output will match the input frame, no matter the coordinate system

        #First we convert rvec to the cartesian frame where the field is defined (rotation axis frame R)
        if self.rotationFlag == "F":
            rvec = np.matmul(self.R, rvec)


        #Assume spherical coordinate system
        #we want positions in r, theta, phi 

        r = np.sqrt(np.power(rvec[0], 2) + np.power(rvec[1], 2) + np.power(rvec[2], 2))
        theta = np.arctan2(np.sqrt(np.power(rvec[0], 2) + np.power(rvec[1], 2)), rvec[2])
        phi = np.arctan2(rvec[1], rvec[0])

        #Analytical form given by Connerney (1993)

        #Will get B as a vector containing Br, Btheta, Bphi

        Br = 0
        Btheta = 0
        Bphi = 0

        frac = self.a/r

        for n in range(1, self.nMax+1):
            frac2 = np.power(frac, (n+2))
            for m in range(0, n+1):
    
                #g[n,m] difference in matrix index to n value
              
                Br += (n+1)*frac2*(self.g[n-1, m]*np.cos(m*phi) + self.h[n-1, m]*np.sin(m*phi))*self.PnmCos(n, m, theta)
                Btheta -= (frac2*(self.g[n-1, m]*np.cos(m*phi) + self.h[n-1, m]*np.sin(m*phi))*self.PnmCosDerivative(n, m, theta))
                Bphi += m*frac2*(self.g[n-1, m]*np.sin(m*phi) - self.h[n-1, m]*np.cos(m*phi))*self.PnmCos(n, m, theta)
        
        Bphi *= 1/np.sin(theta)

        #Now need to recover cartesian mag field components
        Bspherical = np.array([Br, Btheta, Bphi])

        if returnR:
            #In this case, we ignore the rotation flag
            if returnCartesian:
                #Get Bcartesian in R frame
                Bcartesian = self.convertPolarToCartesian(Bspherical, theta, phi)
                #Already in R frame
                return Bcartesian
            else:
                return Bspherical
        else:
            #In this case, the rotation flag is respected
            #Get Bcartesian in R frame
            Bcartesian = self.convertPolarToCartesian(Bspherical, theta, phi)
            #Now convert back to the desired frame
            Bcartesian = np.matmul(self.Rinv, Bcartesian)

            if returnCartesian:
                return np.array(Bcartesian)
            else:
                #Get position vector in original frame (F)
                rvec = np.matmul(self.Rinv, rvec)
                #Extract theta and phi from it 
                thetaOriginal = np.arccos(rvec[2]/np.linalg.norm(rvec))
                phiOriginal = np.arctan2(rvec[1], rvec[0])
                #Calculate Bspherical in original frame 
                Bspherical = self.convertCartesianToPolar(Bcartesian, thetaOriginal, phiOriginal, origin=False)

                return Bspherical


    def getGradB(self, rvec, returnCartesian=True):
        #Input rvec will always be cartesian
        #The gradient of the magnetic field will be computed analytically

        #Need the actual field at each point -- may wish to consider incorporating the field to save
        #computational time.
        #THIS IS THE FIELD IN R FRAME, SPHERICAL POLAR
        Bspherical = self.getField(rvec, returnCartesian=False, returnR=True)

        #First we convert rvec to the cartesian frame where the field is defined (rotation axis frame R)
        rvec = np.matmul(self.R, rvec)


        #Assume spherical coordinate system
        #we want positions in r, theta, phi 
        rvec = self.convertCartesianToPolar(rvec, origin = True)
        r = rvec[0]
        theta = rvec[1]
        phi = rvec[2]


        #Nine individual gradients must be computed
        dBrdr = 0
        dBrdtheta = 0
        dBrdphi = 0
        dBthetadr = 0
        dBthetadtheta = 0
        dBthetadphi = 0
        dBphidr = 0
        dBphidtheta = 0
        dBphidphi = 0

        for n in range(1, self.nMax+1):
            #Acquire some common factors for speed of computation
            common_rFactor = np.power((self.a/r), n+2)
            for m in range(0, n+1):
                C = np.cos(m*phi)
                S = np.sin(m*phi)
                common_phiFactor1 = self.g[n-1, m]*C + self.h[n-1, m]*S
                common_phiFactor2 = self.g[n-1, m]*S - self.h[n-1, m]*C
                common_thetaFactor1 = self.PnmCos(n, m, theta)
                common_thetaFactor2 = self.PnmCosDerivative(n, m, theta)
                common_thetaFactor3 = self.PnmCosDerivative2(n, m, theta)

                

                #Now construct the derivatives
                dBrdr -= (n+1)*(n+2)*self.a**(n+2)*r**(-(n+3))*common_phiFactor1*common_thetaFactor1
                dBrdtheta += (n+1)*common_rFactor*common_phiFactor1*common_thetaFactor2
                dBrdphi += (n+1)*common_rFactor*m*(-self.g[n-1, m]*S + self.h[n-1, m]*C)*common_thetaFactor1
                dBthetadr += (n+2)*self.a**(n+2)*r**(-(n+3))*common_phiFactor1*common_thetaFactor2
                dBthetadtheta -= common_rFactor*common_phiFactor1*common_thetaFactor3
                dBthetadphi -= common_rFactor*m*(-self.g[n-1, m]*S + self.h[n-1, m]*C)*common_thetaFactor2
                dBphidr -= m*(n+2)*self.a**(n+2)*r**(-(n+3))*common_phiFactor2*common_thetaFactor1
                dBphidtheta += m*common_rFactor*common_phiFactor2*(np.sin(theta)*common_thetaFactor2-np.cos(theta)*common_thetaFactor1)/(np.power(np.sin(theta), 2))
                dBphidphi += m*common_rFactor*m*common_phiFactor1*common_thetaFactor1

        dBphidr = dBphidr/np.sin(theta)
        dBphidphi = dBphidphi/np.sin(theta)

        # print(dBrdr)
        # print(dBrdtheta)
        # print(dBrdphi)
        # print(dBthetadr)
        # print(dBthetadtheta)
        # print(dBthetadphi)
        # print(dBphidr)
        # print(dBphidtheta)
        # print(dBphidphi)


        #Now we can assemble the next gradients
        normSpherical = np.linalg.norm(Bspherical)
        dBdr = (Bspherical[0]*dBrdr + Bspherical[1]*dBthetadr + Bspherical[2]*dBphidr)/normSpherical
        dBdtheta = (Bspherical[0]*dBrdtheta + Bspherical[1]*dBthetadtheta + Bspherical[2]*dBphidtheta)/normSpherical
        dBdphi = (Bspherical[0]*dBrdphi + Bspherical[1]*dBthetadphi + Bspherical[2]*dBphidphi)/normSpherical

        gradB = np.array([dBdr, (1/r)*dBdtheta, (1/(r*np.sin(theta)))*dBdphi])

        if not(returnCartesian):
            raise Exception("spherical polar output of gradient not yet possible")
        
        #Get gradient of field in cartesian (R frame)
        gradB = self.convertPolarToCartesian(gradB, theta, phi)
        
        #Now convert to desired frame
        gradB = np.matmul(self.Rinv, gradB)

        return gradB


    def getBCrossGradB(self, rvec):
        B = self.getField(rvec)
        gradB = self.getGradField(rvec)

        return np.cross(B, gradB)
        
    
    def getComponents(a, b):
        #a and b are both vectors
        #Vector components of a parallel to b and perpendicular to b are returned
        para = np.dot(a, b)*b/np.linalg.norm(b)
        perp = a - para
        return para, perp
        

    def getLongitudePlaneB(self, rMax, phi, N, planetaryFilter = True):
        #rMax is the maximum r to plot out to in  units of self.a
        #phi is the particular longitude to plot at in radians
        #planetaryFilter tells teh function whether or not to plot where the planet is
        #N is the  number of points on the side length of the rMax x rMax grid, inclusive of endpoints
        
        B = np.zeros((2*N, N, 3))
        r = np.zeros(np.shape(B))
        
        l = self.a*rMax/(N-1)
        latticeVecAxial = np.array([0, 0, 2*self.a*rMax/(2*N-1)])
        latticeVecRho = np.array([l*np.cos(phi), l*np.sin(phi), 0])
        # startVec = -0.5*(2*N - 1)*latticeVecAxial
        startVec = -(N - 0.5)*latticeVecAxial
        for i in range(0, np.shape(B)[0]):
            for j in range(0, np.shape(B)[1]):
                current_r = startVec + i*latticeVecAxial + j*latticeVecRho
                r[i, j] = current_r
                if planetaryFilter:
                    if np.linalg.norm(current_r) < self.a:
                        B[i, j] = np.zeros((3))
                    else:
                        B[i, j] = self.getField(current_r)
                else:
                    B[i, j] = self.getField(current_r)
        
        return B, r


    def getLongitudePlaneGradB(self, rMax, phi, N, planetaryFilter = True):
        #rMax is the maximum r to plot out to in  units of self.a
        #phi is the particular longitude to plot at in radians
        #planetaryFilter tells teh function whether or not to plot where the planet is
        #N is the  number of points on the side length of the rMax x rMax grid, inclusive of endpoints
        
        gradB = np.zeros((2*N, N, 3))
        r = np.zeros(np.shape(gradB))
        
        l = self.a*rMax/(N-1)
        latticeVecAxial = np.array([0, 0, 2*self.a*rMax/(2*N-1)])
        latticeVecRho = np.array([l*np.cos(phi), l*np.sin(phi), 0])
        # startVec = -0.5*(2*N - 1)*latticeVecAxial
        startVec = -(N - 0.5)*latticeVecAxial
        for i in range(0, np.shape(gradB)[0]):
            for j in range(0, np.shape(gradB)[1]):
                current_r = startVec + i*latticeVecAxial + j*latticeVecRho
                r[i, j] = current_r
                if planetaryFilter:
                    if np.linalg.norm(current_r) < self.a:
                        gradB[i, j] = np.zeros((3))
                    else:
                        gradB[i, j] = self.getGradB(current_r)
                else:
                    gradB[i, j] = self.getGradB(current_r)
        
        return gradB, r


    def getLongitudePlaneDriftDirection(self, rMax, phi, N):
        #First get B
        BPlane, r = self.getLongitudePlaneB(rMax, phi, N)
        
        #Now get gradB
        gradBPlane, r = self.getLongitudePlaneGradB(rMax, phi, N)

        #Now have to combine the two to get the drift direction
        drift = np.cross(BPlane, gradBPlane)
        #Normalise
        drift = drift/np.linalg.norm(drift, axis=-1, keepdims=True)
        
        #nans may emerge so remove - TODO check origin of the nans
        drift = np.nan_to_num(drift)

        return drift, r

    
    def getLongitudePlanesB(self, deltaPhi, rMax, N, planetaryFilter = True):
        #deltaPhi in degrees
        deltaPhi = np.abs(deltaPhi)
        deltaPhi = np.pi/180 * deltaPhi
        noPlanes = int(np.round(2*np.pi/deltaPhi))
        planeBData = np.zeros((noPlanes, 2*N, N, 3))
        planerData = np.zeros(np.shape(planeBData))
        for p in range(0, noPlanes):
            phi = p*deltaPhi
            planeBData[p], planerData[p] = self.getLongitudePlaneB(rMax, phi, N, planetaryFilter)
            
        return planeBData, planerData
        

    def plotLongitudePlanesB(self, deltaPhi, rMax, N, planetaryFilter = True, animate = False):
        #This function should be plotting the magnetic field at multiple longitudes.
        #It isn't
        #You might want to fix that at some point.
        quadrupolePlanes, positionData = self.getLongitudePlanesB(deltaPhi, rMax, N, planetaryFilter=True)
        
        counter = 0

        def getFrame(fig, data, counter):
            fig.imshow(data[counter])
            return
        

        fig = plt.figure(4)
        getFrame(fig, quadrupolePlanes, counter)
        ani = FuncAnimation(fig, getFrame, frames = np.arange(0, 360 + deltaPhi, deltaPhi))

        plt.show()


        return


    def plotDeviationData(self, deltaPhi, rMax, N):
        #This function plots data about how the field deviates from that of a dipole.
        #To do this, we temporarily reset self.nMax

        #deltaPhi in degrees

        #Assume the field has been rotated to F, but it shouldn't be a problem if it's in R

        #First get the dipole data
        true_nMax = np.copy.deepcopy(self.nMax)
        self.nMax = 1

        dipolePlanes, positionData = self.getLongitudePlanesB(deltaPhi, rMax, N, planetaryFilter=True)

        #Now set nMax back to the full field
        self.nMax = np.copy.deepcopy(true_nMax)

        quadrupolePlanes, positionData2 = self.getLongitudePlanesB(deltaPhi, rMax, N, planetaryFilter=True)

        #Now can do maths on the planes to figure out where the important features should be
        
        #Find difference between field with quadrupole and field without
        differencePlanes = quadrupolePlanes - dipolePlanes
        print(differencePlanes)
        #Get magnitude of this quadrupole only field
        differencePlanesMag = np.linalg.norm(differencePlanes, axis = 3)
        #Get magnitude of dipole only field
        dipolePlanesMag = np.linalg.norm(dipolePlanes, axis = 3)
        
        #Find ratio at each point
        ratioPlanes = differencePlanesMag/dipolePlanesMag
        #Replace any nan with 0.0
        ratioPlanes = np.nan_to_num(ratioPlanes) #replace nan with 0.0

        #Find maxiumum on each axis until we have the maximum at each 
        maxRatioQ_D = np.amax(ratioPlanes, axis = 2)
        maxRatioQ_D = np.amax(maxRatioQ_D, axis = 1)

        ax3 = plt.figure().add_subplot()
        xAxis = np.arange(0, 360, deltaPhi)
        ax3.plot(xAxis, maxRatioQ_D)
        ax3.set_xlabel("Longitude ($\degree$)")
        ax3.set_ylabel("Maximum ratio of absolute deviation from dipole against dipole field")
        ax3.set_title("Field Deviation")


    def plotDeviationColourMapLongitudePlane(self, phi, rMax, N, vectorStep = 10, planetaryFilter = True, plot=True):
        #This plots the deviation of the complete field from a dipole for a single plane of constant longitude
        #vecPerRadius is the number of vectors to be plotted per planet radius along an axis

        #First some input checks on phi:
        if isinstance(phi, int):
            phiArray = np.array([phi])
        elif isinstance(phi, list):
            phiArray = np.array(phi)
        else:
            phiArray = copy.deepcopy(phi)

        
        
        diagnostic = np.zeros((np.shape(phiArray)[0], 2*N, N))
        vectors = np.zeros((np.shape(phiArray)[0], 2*N, N, 2))
        vecPos = np.zeros(np.shape(vectors))
        counter = 0

        phiArrayDeg = copy.deepcopy(phiArray)
        phiArray = np.deg2rad(phiArray)

        for phi in phiArray:
            #Get complete B field for a particular plane
            completePlane, positionData = self.getLongitudePlaneB(rMax, phi, N, planetaryFilter=planetaryFilter)
            #Get dipole only B field for a particular plane
            true_nMax = copy.deepcopy(self.nMax)
            self.nMax = 1
            dipolePlane, positionData = self.getLongitudePlaneB(rMax, phi, N, planetaryFilter=planetaryFilter)
            self.nMax = copy.deepcopy(true_nMax)

            #Now find the difference between the fields
            quadrupolePlane = completePlane - dipolePlane
            #Calculate the diagnostic
            diagnostic[counter] = np.nan_to_num(np.linalg.norm(quadrupolePlane, axis = 2)/np.linalg.norm(dipolePlane, axis = 2))
            #Now need to calculate the vectors to be plotted

            a = np.array([np.cos(phi), np.sin(phi), 0])
            b = np.array([0, 0, 1])
            planeVec = np.array([np.dot(completePlane, a), np.dot(completePlane, b)])
            #The dimension layout of the above array is disgusting
            vectors[counter, :, :, 0] = planeVec[0]
            vectors[counter, :, :, 1] = planeVec[1]

            vecPos[counter, :, :, 0] = np.sqrt(np.power(positionData[:, :, 0], 2) + np.power(positionData[:, :, 1], 2))
            vecPos[counter, :, :, 1] = positionData[:, :, 2]
            
            counter+=1
        
        vectors = vectors/np.linalg.norm(vectors, axis=3, keepdims=True)
        vecPos = vecPos/self.a

        #Vectors are usually plotted at a lesser density than the colour map
        #every vectorStep'th value is plotted

        vectors2 = np.zeros((np.shape(vectors)[0], int(np.ceil(np.shape(vectors)[1]/(vectorStep))), int(np.ceil(np.shape(vectors)[2]/vectorStep)), np.shape(vectors)[3]))
        vecPos2 = np.zeros((np.shape(vectors)[0], int(np.ceil(np.shape(vectors)[1]/(vectorStep))), int(np.ceil(np.shape(vectors)[2]/vectorStep)), np.shape(vectors)[3]))

        for i in range(0, np.shape(vectors2)[1]):
            for j in range(0, np.shape(vectors2)[2]):
                vecPos2[:, i, j, :] = vecPos[:, vectorStep*i, vectorStep*j, :]
                vectors2[:, i, j, :] = vectors[:, vectorStep*i, vectorStep*j, :]

        #Finally, need to create the curves showing constant L shell
        #These are independent of phi
        LShells = np.arange(0, rMax, 2)
        L_theta = np.arange(0, np.pi, 0.01)
        L_x = np.zeros((np.shape(LShells)[0], np.shape(L_theta)[0]))
        L_y = np.zeros(np.shape(L_x))
        for i in range(0, np.shape(LShells)[0]):
            # L_r[i, :] = LShells[i]*np.power(np.sin(L_theta), 2)
            L_x[i, :] = LShells[i]*np.power(np.sin(L_theta), 3)
            L_y[i, :] = LShells[i]*np.power(np.sin(L_theta), 2)*np.cos(L_theta)


        #Now on to the plotting
        phiArray = np.rad2deg(phiArray)
       

        if plot:
            noFigs = np.shape(phiArray)[0]
            maxCols = 6
            figCols = int(min([maxCols, noFigs]))
            figRows = int(np.ceil(noFigs/figCols))
            fig, axs = plt.subplots(nrows=figRows, ncols=figCols, squeeze=False)
            rhoAxis = np.linspace(0, rMax, N)
            zAxis = np.linspace(-rMax, rMax, 2*N)

            flatDiagnostic = diagnostic.flatten()
            flatDiagnostic = np.sort(flatDiagnostic)
            found = False
            counter = 0
            while found == False:
                if flatDiagnostic[counter] != 0.0:
                    found = True
                    colourMin = flatDiagnostic[counter]
                else:
                    counter+=1

            colourMax = np.max(flatDiagnostic)

            print(colourMin)
            print(colourMax)
            
            counter = 0
            for i in range(0, figRows):
                for j in range(0, figCols):
                    if counter < noFigs:
                        # obj = axs[i].pcolormesh(rhoAxis, zAxis, diagnostic, cmap = "plasma", norm=LogNorm())
                        obj = axs[i, j].pcolormesh(rhoAxis, zAxis, diagnostic[counter, :, :], cmap = "plasma", norm=LogNorm(), vmax=colourMax, vmin=colourMin)
                        axs[i, j].quiver(vecPos2[counter, :, :, 0], vecPos2[counter, :, :, 1], vectors2[counter, :, :, 0], vectors2[counter, :, :, 1])
                        for k in range(0, np.shape(L_x)[0]):
                            axs[i, j].plot(L_x[k,:], L_y[k,:], color = "black")
                        axs[i, j].set_xlabel("rho")
                        axs[i, j].set_ylabel("z")
                        titleString = "Phi = " + str(phiArrayDeg[counter]) + "$\degree$"
                        axs[i, j].set_title(titleString)
                        axs[i, j].set_aspect("equal") 
                        counter += 1
            
            fig.colorbar(obj, ax = axs.ravel().tolist())
        
        return
        

    # def getLShellB(self, L, NTheta = 100, NPhi = 361):
    #     #Returns an array of the magnetic field at a particular L-shell 
    #     #deltaTheta in degrees
    #     #deltaPhi in degrees

    #     #First get the arrays set up
    #     #Calculate cut-off angle
    #     cutOffTheta = np.arcsin(np.sqrt(1/L))
    #     thetaArray = np.linspace(cutOffTheta, np.pi - cutOffTheta, NTheta)
    #     deltaPhi = 2*np.pi/(NPhi-1)

    #     #There are Ny rows, Nx columns, and each point has 3 things, r theta and phi
    #     rF_spherical = np.zeros((NTheta, NPhi, 3))
    #     rF_cartesian = np.zeros(np.shape(rF_spherical))
    #     bF = np.zeros(np.shape(rF_cartesian))
    #     for i in range(0, NTheta):
    #         for j in range(0, NPhi):
    #             #set theta
    #             rF_spherical[i, j, 1] = np.copy(thetaArray[i])
    #             #set phi
    #             rF_spherical[i, j, 2] = deltaPhi * j
    #             #set r
    #             rF_spherical[i, j, 0] = self.a*L*np.power(np.sin(rF_spherical[i, j, 1]), 2)

    #              #Now have the coordinates in field-aligned spherical polars
    #             #Need to convert to field aligned cartesian

    #             theta = copy.deepcopy(rF_spherical[i, j, 1])
    #             phi = copy.deepcopy(rF_spherical[i, j, 2])
    #             rF_spherical[i, j, 1] = 0
    #             rF_spherical[i, j, 2] = 0
    #             rF_cartesian[i, j, :] = self.convertPolarToCartesian(rF_spherical[i, j, :], theta, phi)

    #             #Now get B field at each point
    #             bF[i, j, :] = self.getField(rF_cartesian[i, j, :])

    #             rF_spherical[i, j, 1] = copy.copy(theta)
    #             rF_spherical[i, j, 2] = copy.copy(phi)
 
               
    #     return rF_spherical, rF_cartesian, bF

    def getLShellB(self, L, deltaTheta = 1, deltaPhi = 1):
        #Calculate cut-off angle
        cutOffTheta = abs(np.arcsin(np.sqrt(1/L)))
        thetaArray = np.arange(0, np.pi + deltaTheta*np.pi/180, deltaTheta*np.pi/180)
        phiArray = np.arange(0, 2*np.pi + deltaPhi*np.pi/180, deltaPhi*np.pi/180)

        positions = np.zeros((np.shape(thetaArray)[0], np.shape(phiArray)[0], 3))
        B = np.zeros(np.shape(positions))

        #calculate positions
        for i in range(0, len(thetaArray)):
            r = self.a*L*np.sin(thetaArray[i])**2
            for j in range(0, len(phiArray)):
                positions[i, j, :] = self.convertPolarToCartesian(np.array([r, 0, 0]), thetaArray[i], phiArray[j])
                if (thetaArray[i] >= cutOffTheta) and (thetaArray[i] <= (np.pi - cutOffTheta)):
                    B[i, j, :] = self.getField(positions[i, j, :])
                else:
                    B[i, j, :] = 0
        
        return positions, B

    def getLShellGradB(self, L, deltaTheta = 1, deltaPhi = 1):
        #Calculate cut-off angle
        cutOffTheta = np.arcsin(np.sqrt(1/L))
        thetaArray = np.arange(0, np.pi + deltaTheta*np.pi/180, deltaTheta*np.pi/180)
        phiArray = np.arange(0, 2*np.pi + deltaPhi*np.pi/180, deltaPhi*np.pi/180)

        positions = np.zeros((np.shape(thetaArray)[0], np.shape(phiArray)[0], 3))
        gradB = np.zeros(np.shape(positions))

        #calculate positions
        for i in range(0, len(thetaArray)):
            r = self.a*L*np.sin(thetaArray[i])**2
            for j in range(0, len(phiArray)):
                positions[i, j, :] = self.convertPolarToCartesian(np.array([r, 0, 0]), thetaArray[i], phiArray[j])
                if (thetaArray[i] >= cutOffTheta) and (thetaArray[i] <= (np.pi - cutOffTheta)):
                    gradB[i, j, :] = self.getGradB(positions[i, j, :])
                else:
                    gradB[i, j, :] = 0
        
        return positions, gradB

    def getLShellDriftDirection(self,  L, deltaTheta = 1, deltaPhi = 1):
        positions, B = self.getLShellB(L, deltaTheta=deltaTheta, deltaPhi=deltaPhi)
        positions, gradB = self.getLShellGradB(L, deltaTheta=deltaTheta, deltaPhi=deltaPhi)

    

        #Compute the drift direction
        drift = np.cross(B, gradB)
        #Normalise
        drift = drift/np.linalg.norm(drift, axis=-1, keepdims=True)
        

        return drift, positions


    
    def plotDeviationColourMapLShell(self, L, NTheta = 100, NPhi = 361, vectorStep = 10, plot=True):
        #L is the L shell to plot at
        #deltaTheta is the interval in degrees between points
        #deltaPhi is the interval in degrees between points in the phi direction
        

        #returns a plot

        if L < 1:
            raise(Exception("Warning!  You're attempting to plot at an L-shell inside the planet!"))

        #Get L-shell data for complete field
      
        rF_spherical, rF_cartesian, LShellData_complete = self.getLShellB(L, NTheta = NTheta, NPhi = NPhi)
        
        #change nMax to get dipole field only
        true_nMax = np.copy(self.nMax)
        self.nMax = 1
        rF_spherical, rF_cartesian, LShellData_dipole = self.getLShellB(L, NTheta = NTheta, NPhi = NPhi)
        #change nMax back so it's the full field
        self.nMax = np.copy(true_nMax)

        LShellData_quadrupole = LShellData_complete - LShellData_dipole
        ratio = np.linalg.norm(LShellData_quadrupole, axis = 2)/np.linalg.norm(LShellData_dipole, axis = 2)

        #Now need to compute data for vectors
        #First get an array of deviation angles
        vectors = np.zeros((NTheta, NPhi, 2))
        vectors[:,:,1] = 1
        print(vectors[0, 0, :])
        #vectors is now a grid of (0, 1) vectors pointing up
        for i in range(0, NTheta):
            for j in range(0, NPhi):
                gamma = np.arccos(((np.dot(LShellData_complete[i, j, :], LShellData_dipole[i, j, :]))/(np.linalg.norm(LShellData_complete[i, j, :])*np.linalg.norm(LShellData_dipole[i, j, :]))))
                vectors[i, j, :] = self.rotate2D(vectors[i, j, :], gamma)

        vecPos = np.zeros(np.shape(vectors))
        vecPos[:, :, 0] = rF_spherical[:, :, 1]
        vecPos[:, :, 1] = rF_spherical[:, :, 2]
        vecPos = vecPos*180/np.pi

        if plot:
            vectors2 = np.zeros((int(np.ceil(np.shape(vectors)[0]/vectorStep)), int(np.ceil(np.shape(vectors)[1]/vectorStep)), np.shape(vectors)[2]))
            vecPos2 = np.zeros(np.shape(vectors2))
            for i in range(0, np.shape(vectors2)[0]):
                for j in range(0, np.shape(vectors2)[1]):
                    vecPos2[i, j, :] = vecPos[vectorStep*i, vectorStep*j, :]
                    vectors2[i, j, :] = vectors[vectorStep*i, vectorStep*j, :]
                
                
                
            
            
            thetaAxis = rF_spherical[:, 0, 1]*180/np.pi
            phiAxis = rF_spherical[0, :, 2]*180/np.pi
            # phiAxis = np.arange(0, 360 + deltaPhi, deltaPhi)
            
            ax6 = plt.figure().add_subplot()
            # imshowObject = ax6.imshow(np.transpose(ratio), cmap = "plasma", norm=LogNorm(vmin=0.1, vmax=1))
            # imshowObject = ax6.imshow(np.transpose(ratio), cmap = "plasma", norm=LogNorm())
            obj = ax6.pcolormesh(phiAxis, thetaAxis, ratio, cmap = "plasma", norm=LogNorm())
            ax6.quiver(vecPos2[:, :, 1], vecPos2[:, :, 0], vectors2[:, :, 0], vectors2[:, :, 1])
            ax6.set_ylim(ax6.get_ylim()[::-1])
            plt.colorbar(obj)
            ax6.set_xlabel("Phi ($\degree$)")
            ax6.set_ylabel("Theta ($\degree$)")
            titleString = "Diagnostic on L = " + str(L)
            ax6.set_title(titleString)
            ax6.set_aspect("equal")

        return


    def plotDeviationColourMapLShell2(self, L, NTheta = 100, NPhi = 361, vectorStep = 10, plot=True):
        #UNFINISHED:  purpose is to plot multiple L shell graphs as subplots.  THe plots change size, so you'll have
        #to find a way to deal with that

        #L is the L shell to plot at
        #deltaTheta is the interval in degrees between points
        #deltaPhi is the interval in degrees between points in the phi direction

        #First some input checks on phi:
        if isinstance(L, int):
            LArray = np.array([L])
        elif isinstance(L, list):
            LArray = np.array(L)
        else:
            LArray = copy.deepcopy(L)

        diagnostic = np.zeros((np.shape(LArray)[0], NTheta, NPhi))
        vectors = np.zeros((np.shape(LArray)[0], NTheta, NPhi, 2))
        vectors[:, :, :, 1] = 1
        vecPos = np.zeros(np.shape(vectors))

        counter = 0
        for L in LArray:
            if L < 1:
                raise(Exception("Warning!  You're attempting to plot at an L-shell inside the planet!"))

            #Get L-shell data for complete field
        
            rF_spherical, rF_cartesian, LShellData_complete = self.getLShellB(L, NTheta = NTheta, NPhi = NPhi)
            
            #change nMax to get dipole field only
            true_nMax = np.copy(self.nMax)
            self.nMax = 1
            rF_spherical, rF_cartesian, LShellData_dipole = self.getLShellB(L, NTheta = NTheta, NPhi = NPhi)
            #change nMax back so it's the full field
            self.nMax = np.copy(true_nMax)

            LShellData_quadrupole = LShellData_complete - LShellData_dipole
            diagnostic[counter] = np.linalg.norm(LShellData_quadrupole, axis = 2)/np.linalg.norm(LShellData_dipole, axis = 2)

            #Now need to compute data for vectors
            #First get an array of deviation angles
            

            #vectors is now a grid of (0, 1) vectors pointing up
            for i in range(0, NTheta):
                for j in range(0, NPhi):
                    gamma = np.arccos(((np.dot(LShellData_complete[i, j, :], LShellData_dipole[i, j, :]))/(np.linalg.norm(LShellData_complete[i, j, :])*np.linalg.norm(LShellData_dipole[i, j, :]))))
                    vectors[counter, i, j, :] = self.rotate2D(vectors[i, j, :], gamma)

            vecPos[counter, :, :, 0] = rF_spherical[:, :, 1]
            vecPos[counter, :, :, 1] = rF_spherical[:, :, 2]
            
            counter+=1

        vecPos = vecPos*180/np.pi

        if plot:
            vectors2 = np.zeros((np.shape(LArray)[0], int(np.ceil(np.shape(vectors)[0]/vectorStep)), int(np.ceil(np.shape(vectors)[1]/vectorStep)), np.shape(vectors)[2]))
            vecPos2 = np.zeros(np.shape(vectors2))
            for i in range(0, np.shape(vectors2)[1]):
                for j in range(0, np.shape(vectors2)[2]):
                    vecPos2[:, i, j, :] = vecPos[:, vectorStep*i, vectorStep*j, :]
                    vectors2[:, i, j, :] = vectors[:, vectorStep*i, vectorStep*j, :]
                
                
            noFigs = np.shape(LArray)[0]
            maxCols = 4
            figCols = int(min([maxCols, noFigs]))
            figRows = int(np.ceil(noFigs/figCols))
            fig, axs = plt.subplots(nrows=figRows, ncols=figCols, squeeze=False)    
            
            
            thetaAxis = rF_spherical[:, 0, 1]*180/np.pi
            phiAxis = rF_spherical[0, :, 2]*180/np.pi
            # phiAxis = np.arange(0, 360 + deltaPhi, deltaPhi)
            counter = 0
            for i in range(0, figRows):
                 for j in range(0, figCols):
                     if counter < noFigs:
                        # imshowObject = ax6.imshow(np.transpose(ratio), cmap = "plasma", norm=LogNorm(vmin=0.1, vmax=1))
                        # imshowObject = ax6.imshow(np.transpose(ratio), cmap = "plasma", norm=LogNorm())
                        obj = axs[i, j].pcolormesh(phiAxis, thetaAxis, diagnostic, cmap = "plasma", norm=LogNorm())
                        axs[i, j].quiver(vecPos2[:, :, 1], vecPos2[:, :, 0], vectors2[:, :, 0], vectors2[:, :, 1])
                        axs[i, j].set_ylim(axs[i, j].get_ylim()[::-1])
                        axs[i, j].set_xlabel("Phi ($\degree$)")
                        axs[i, j].set_ylabel("Theta ($\degree$)")
                        titleString = "Diagnostic on L = " + str(L)
                        axs[i, j].set_title(titleString)
                        axs[i, j].set_aspect("equal")

            fig.colorbar(obj, ax = axs.ravel().tolist())


            #  counter = 0
            # for i in range(0, figRows):
            #     for j in range(0, figCols):
            #         if counter < noFigs:
            #             # obj = axs[i].pcolormesh(rhoAxis, zAxis, diagnostic, cmap = "plasma", norm=LogNorm())
            #             obj = axs[i, j].pcolormesh(rhoAxis, zAxis, diagnostic[counter, :, :], cmap = "plasma", norm=LogNorm(), vmax=colourMax, vmin=colourMin)
            #             axs[i, j].quiver(vecPos2[counter, :, :, 0], vecPos2[counter, :, :, 1], vectors2[counter, :, :, 0], vectors2[counter, :, :, 1])
            #             for k in range(0, np.shape(L_x)[0]):
            #                 axs[i, j].plot(L_x[k,:], L_y[k,:], color = "black")
            #             axs[i, j].set_xlabel("rho")
            #             axs[i, j].set_ylabel("z")
            #             titleString = "Phi = " + str(phiArrayDeg[counter]) + "$\degree$"
            #             axs[i, j].set_title(titleString)
            #             axs[i, j].set_aspect("equal") 
            #             counter += 1
            
        



        return


    def plotDriftDirectionLongitudePlane(self, phi, rMax, N, plot=True, vectorStep = 10):
        #phi in degrees


        #Check field orientation
        if self.rotationFlag == "R":
            raise(Exception("The field must be rotated first!"))
        #First some input checks on phi:
        if isinstance(phi, int):
            phiArray = np.array([phi])
        elif isinstance(phi, list):
            phiArray = np.array(phi)
        else:
            phiArray = copy.deepcopy(phi)

        phiArrayDeg = copy.deepcopy(phiArray)
        phiArray = np.deg2rad(phiArray)

        drift = np.zeros((np.shape(phiArray)[0], 2*N, N, 3)) #drift direction vectors
        positions = np.zeros(np.shape(drift)) #position vectors cartesian F frame
        normalVec = np.zeros(np.shape(drift))
        counter = 0
        for phi in phiArray:
            #Obtain plane of drift vectors
            drift[counter, :, :, :], positions[counter, :, :, :] = self.getLongitudePlaneDriftDirection(rMax, phi, N)
            #Obtain the cartesian components of each plane's normal vector
            normalVec[counter, :, :, 0] = np.sin(phi)
            normalVec[counter, :, :, 1] = -np.cos(phi)
            counter += 1


        #We now have the drift direction vectors in cartesian, oriented with z parallel to dipole axis.
        #Drift vectors are already normalised
        #Also have plane normal vectors.
        
        #Compute normal component for colour map
        normalComponent = np.multiply(drift, normalVec)
        normalComponent = np.sum(normalComponent, axis = -1)
       

        #Now need to compute the vector components within the plane
        planeVecs = np.zeros((np.shape(phiArray)[0], 2*N, N, 2))
        vecPos = np.zeros(np.shape(planeVecs))
        vecPos[:, :, :, 1] = positions[:, :, :, 2]/self.a
        vecPos[:, :, :, 0] = np.sqrt(positions[:, :, :, 0]**2 + positions[:, :, :, 1]**2)/self.a
        counter = 0
        for phi in phiArray:
            planeVecs[counter, :, :, 0] = np.cos(phi)*drift[counter, :, :, 0] + np.sin(phi)*drift[counter, :, :, 1]
            counter += 1
        planeVecs[:, :, :, 1] = drift[:, :, :, 2]

        #Now reduce density of vectors
        #Vectors are usually plotted at a lesser density than the colour map
        #every vectorStep'th value is plotted

        planeVecs2 = np.zeros((np.shape(planeVecs)[0], int(np.ceil(np.shape(planeVecs)[1]/(vectorStep))), int(np.ceil(np.shape(planeVecs)[2]/vectorStep)), np.shape(planeVecs)[3]))
        vecPos2 = np.zeros((np.shape(planeVecs)[0], int(np.ceil(np.shape(planeVecs)[1]/(vectorStep))), int(np.ceil(np.shape(planeVecs)[2]/vectorStep)), np.shape(planeVecs)[3]))

        for i in range(0, np.shape(planeVecs2)[1]):
            for j in range(0, np.shape(planeVecs2)[2]):
                vecPos2[:, i, j, :] = vecPos[:, vectorStep*i, vectorStep*j, :]
                planeVecs2[:, i, j, :] = planeVecs[:, vectorStep*i, vectorStep*j, :]



        #Now plot

        if plot:
            noFigs = np.shape(phiArray)[0]
            maxCols = 6
            figCols = int(min([maxCols, noFigs]))
            figRows = int(np.ceil(noFigs/figCols))
            fig, axs = plt.subplots(nrows=figRows, ncols=figCols, squeeze=False)
            rhoAxis = np.linspace(0, rMax, N)
            zAxis = np.linspace(-rMax, rMax, 2*N)

            flatNormal = normalComponent.flatten()
            flatNormal = np.sort(flatNormal)
            found = False
            counter = 0
            while found == False:
                if flatNormal[counter] != 0.0:
                    found = True
                    colourMin = flatNormal[counter]
                else:
                    counter+=1

            colourMax = np.max(flatNormal)

            #diagnostic pixel
            # normalComponent[0, 0, 0] = -1
            # planeVecs2[0, 0, 0, 1] = 10
            
            counter = 0
            for i in range(0, figRows):
                for j in range(0, figCols):
                    if counter < noFigs:
                        # obj = axs[i].pcolormesh(rhoAxis, zAxis, diagnostic, cmap = "plasma", norm=LogNorm())
                        obj = axs[i, j].pcolormesh(rhoAxis, zAxis, normalComponent[counter, :, :], cmap = "plasma", vmin = -1, vmax = 1)#, norm=LogNorm(vmin=colourMin, vmax=colourMax))
                        axs[i, j].quiver(vecPos2[counter, :, :, 0], vecPos2[counter, :, :, 1], planeVecs2[counter, :, :, 0], planeVecs2[counter, :, :, 1])
                        
                        axs[i, j].set_xlabel("rho")
                        axs[i, j].set_ylabel("z")
                        titleString = "Phi = " + str(phiArrayDeg[counter]) + "$\degree$"
                        axs[i, j].set_title(titleString)
                        axs[i, j].set_aspect("equal") 
                        counter += 1
            
            fig.colorbar(obj, ax = axs.ravel().tolist())
            fig.suptitle("Drift direction vector on slices of constant longitude")




        return


    def plotDriftDirectionLShell(self, LArray, deltaTheta = 1, deltaPhi = 1, vectorStep = 10):
        #Drift will be plotted outward from the equator
        if isinstance(LArray, int):
            LArray = np.array([LArray])
        elif isinstance(LArray, list):
            LArray = np.array(LArray)
        else:
            LArray = copy.deepcopy(LArray)

        thetaArray = np.arange(0, np.pi + deltaTheta*np.pi/180, deltaTheta*np.pi/180)
        phiArray = np.arange(0, 2*np.pi + deltaPhi*np.pi/180, deltaPhi*np.pi/180)
        driftPlanes = np.zeros((np.shape(LArray)[0], np.shape(thetaArray)[0], np.shape(phiArray)[0], 3))
        driftPositions = np.zeros(np.shape(driftPlanes))

        counter = 0
        for L in LArray:
            driftPlanes[counter, :, :, :], driftPositions[counter, :, :, :] = self.getLShellDriftDirection(L, deltaTheta=deltaTheta, deltaPhi=deltaPhi)
            counter += 1
        
        #Compute the normal vectors at each point (not normalised)
        xxHat = copy.deepcopy(driftPositions)
        xxHat[:, :, :, 1] = 0
        xxHat[:, :, :, 2] = 0
        yyHat = copy.deepcopy(driftPositions)
        yyHat[:, :, :, 0] = 0
        yyHat[:, :, :, 2] = 0
        zzHat = copy.deepcopy(driftPositions)
        zzHat[:, :, :, 0] = 0
        zzHat[:, :, :, 1] = 0

        secondFactor = (driftPositions[:, :, :, 2]**2)/(driftPositions[:, :, :, 0]**2 + driftPositions[:, :, :, 1]**2)
        secondFactor = np.array([secondFactor, secondFactor, secondFactor])

        secondFactor = np.transpose(secondFactor, (1, 2, 3, 0))
        normalVecs = driftPositions[:,:,:,:] -2*secondFactor*(xxHat + yyHat) + 2*zzHat  
        #Now normalise them
        normalVecs = normalVecs/np.linalg.norm(normalVecs, axis=-1, keepdims=True)
        #Extract drift components
        normalDrifts = np.multiply(driftPlanes, normalVecs)
        normalDrifts = np.sum(normalDrifts, axis=-1)
        normalDrifts2 = np.array([normalDrifts, normalDrifts, normalDrifts])
        normalDrifts2 = np.transpose(normalDrifts2, (1, 2, 3, 0))
        inPlaneDrifts = driftPlanes - normalDrifts2*normalVecs
        #Now need to convert the 3D in-plane drifts into their 2D equivalents 
        #Define phi hat direction
        phiFactor = 1/(np.sqrt(driftPositions[:, :, :, 0]**2 + driftPositions[:, :, :, 1]**2))
        phiFactor = np.array([phiFactor, phiFactor, phiFactor])
        phiHat = (np.array([-driftPositions[:, :, :, 1], driftPositions[:, :, :, 0], np.zeros(np.shape(driftPositions[:, :, :, 0]))]))/phiFactor
 
        
        phiHat = np.transpose(phiHat, (1, 2, 3, 0))
        #Define the sigma hat direction to be in the plane, perpendicular to the normal vec and the phi hat direction, "going downwards in same-ish direction as theta hat"

        sigmaHat = np.cross(phiHat, normalVecs)
        v = -np.sum(np.multiply(inPlaneDrifts, sigmaHat), axis = -1) #-alpha
        u  = np.sum(np.multiply(inPlaneDrifts, phiHat), axis=-1) #beta
        vectorData = np.zeros((np.shape(LArray)[0], np.shape(thetaArray)[0], np.shape(phiArray)[0], 2))
        vectorData[:, :, :, 0] = u
        vectorData[:, :, :, 1] = v      

        colourMin = -1
        colourMax = 1
        print(colourMin, colourMax)

        phiArray = np.rad2deg(phiArray)
        thetaArray = np.rad2deg(thetaArray)

        #Now reduce density of vectors
        #Vectors are usually plotted at a lesser density than the colour map
        #every vectorStep'th value is plotted

        vectorData2 = np.zeros((np.shape(vectorData)[0], int(np.ceil(np.shape(vectorData)[1]/(vectorStep))), int(np.ceil(np.shape(vectorData)[2]/vectorStep)), np.shape(vectorData)[3]))
        vec_x = np.zeros((int(np.ceil(np.shape(phiArray)[0]/(vectorStep)))))
        vec_y = np.zeros((int(np.ceil(np.shape(thetaArray)[0]/(vectorStep)))))

        for i in range(0, np.shape(vectorData2)[1]):
            for j in range(0, np.shape(vectorData2)[2]):
                vectorData2[:, i, j, :] = vectorData[:, vectorStep*i, vectorStep*j, :]
        
        for i in range(0, np.shape(vec_x)[0]):
            vec_x[i] = phiArray[vectorStep*i]
        
        for i in range(0, np.shape(vec_y)[0]):
            vec_y[i] = thetaArray[vectorStep*i]

        print(np.shape(vec_x))
        print(np.shape(vec_y))
        print(np.shape(vectorData))
        print(np.shape(vectorData2))
        superTitle = "Drift direction vector on planes of constant L-shell"
        self.generalLShellPlot(LArray, vec_x = vec_x, vec_y=vec_y, vectorData = vectorData2, colour_x=phiArray, colour_y=thetaArray, colourData = normalDrifts, colourMin=colourMin, colourMax=colourMax, title=superTitle)

        return


    def generalLShellPlot(self, LArray, vec_x = 0, vec_y = 0, vectorData = 0, colour_x = 0, colour_y = 0, colourData = 0, colourMin = 0.1, colourMax = 1, title=""):
        #All vectors in cartesian, F frame
        #Vector positions already given relative to the plotting surface
        #Need to extract phi axis and theta axis from 
        print("Warning:  Normalisation needs more options.")
        
        noFigs = np.shape(LArray)[0]
        maxCols = 2
        figCols = int(min([maxCols, noFigs]))
        figRows = int(np.ceil(noFigs/figCols))
        fig, axs = plt.subplots(nrows=figRows, ncols=figCols, squeeze=False)  

        # #Diagnostic pixel
        # colourData[0, 0, 0] = 1
        # vectorData[0, 0, 0, 0] = 10
        # vectorData[0, 0, 0, 1] = -10
            
        counter = 0
        for i in range(0, figRows):
            for j in range(0, figCols):
                if counter < noFigs:
                

                    obj = axs[i, j].pcolormesh(colour_x, colour_y, colourData[counter, :, :], cmap = "plasma", vmin=colourMin, vmax=colourMax)#, norm=LogNorm(), vmax=colourMax, vmin=colourMin)
                    axs[i, j].quiver(vec_x, vec_y, vectorData[counter, :, :, 0], vectorData[counter, :, :, 1])
                    axs[i, j].set_ylim(axs[i, j].get_ylim()[::-1])
                    axs[i, j].set_xlabel("Phi ($\degree$)")
                    axs[i, j].set_ylabel("Theta ($\degree$)")
                    titleString = "L = " + str(LArray[counter])
                    axs[i, j].set_title(titleString)
                    axs[i, j].set_aspect("equal")
                    counter+=1

        fig.colorbar(obj, ax = axs.ravel().tolist())
        fig.suptitle(title)

        return



class OTDField(Field):
    def __init__(self, m):
        self.m = m
        self.naturalUnits = False

    def getField(self, rvec):
        r = np.linalg.norm(rvec)
        return -(1/np.power(r, 3))*self.m + (3/np.power(r, 5))*np.dot(self.m, rvec)*rvec


class UniformField(Field):
    def __init__(self, fieldVec):
        self.B = fieldVec
        self.naturalUnits = False
    
    def getField(self, r):
        return self.B
    
    def setNaturalUnits(self, activate, charge, restMass, period):
        #This switches the units to being natural units
        if activate:
            self.B = (charge*period/restMass)*self.B
  
            self.naturalUnits = True
        else:
            self.B = (restMass/(charge*period))*self.B
           
            self.naturalUnits = False
        return



class UranusField(SHField):
    def __init__(self, dipoleOnly=False):
        Ru = 25600000           #radius of Uranus in metres

        #Create field of Uranus
        g = np.array([[11278, 10928, 0], [-9648, -12284, 1453]]) #these are in nanoteslas
        h = np.array([[0, -16049, 0], [0, 6405, 4220]])
        g = g/1000000000 #now in teslas
        h = h/1000000000

        super(UranusField, self).__init__(Ru, g, h, 0, 0, dipoleOnly=dipoleOnly)
        
        return

class UranusFieldOld(SHField):
    def __init__(self, dipoleOnly=False):
        Ru = 25600000           #radius of Uranus in metres

        #Create field of Uranus
        g = np.array([[11893, 11579, 0], [-6030, -12587, 196]]) #these are in nanoteslas
        h = np.array([[0, -15684, 0], [0, 6116, 4759]])
        g = g/1000000000 #now in teslas
        h = h/1000000000

        super(UranusFieldOld, self).__init__(Ru, g, h, 0, 0, dipoleOnly=dipoleOnly)
        
        return

class NeptuneField(SHField):
    pass