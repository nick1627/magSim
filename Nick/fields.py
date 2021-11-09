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
from numpy.lib.arraysetops import isin
import pyshtools as sh
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LogNorm
import copy

class Field:
    def __init__(self):
        self.B = 0
        
        #     #basic field superclass to deal with coordinate system conversions
        #     #we assume that the field wants to be in cartesian at the base level
        return

    def convertCartesianToPolar(self, rvec):
        #converts vectors from the origin only
        
        r = np.sqrt(rvec[0]**2 + rvec[1]**2 + rvec[2]**2)
        theta = np.arctan2(np.sqrt(rvec[0]**2 + rvec[1]**2), rvec[2])
        phi = np.arctan2(rvec[1], rvec[0])

        return np.array([r, theta, phi])


    def convertPolarToCartesian(self, rvec, theta, phi):
        #rvec is the components of the vector to be converted in order (a, b, c) for a*rhat + b*theta hat + c*phihat
        #For vectors from the origin, rvec should read (r, 0, 0), but theta and phi must be given too

        T = np.array([[np.sin(theta)*np.cos(phi),   np.cos(theta)*np.cos(phi),  -np.sin(phi)], 
                        [np.sin(theta)*np.sin(phi), np.cos(theta)*np.sin(phi),  np.cos(phi)],
                        [np.cos(theta),             -np.sin(theta),             0]])

        return np.matmul(T, rvec)


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

        self.rotationFlag = "R"
        self.R = np.identity(3)
        self.Rinv = np.identity(3)

    def rotate(self, rotationKey):
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


    
    def getField(self, rvec):
        #Input rvec will always be cartesian

        #First we convert rvec to the cartesian frame where the field is defined (rotation axis frame R)
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

        Bcartesian = self.convertPolarToCartesian(Bspherical, theta, phi)

        #Now need to convert back to the desired frame
        Bcartesian = np.matmul(self.Rinv, Bcartesian)

        return np.array(Bcartesian)

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
        ax3.set_xlabel("Longitude (º)")
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

        
        
        diagnostic = np.zeros((np.shape(phi)[0], 2*N, N))
        vectors = np.zeros((np.shape(phi)[0], 2*N, N, 2))
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

        # rowsToKeep = np.arange(0, np.shape(vectors)[1], vectorStep)
        # colsToKeep = np.arange(0, np.shape(vectors)[2], vectorStep)
        # colsToDelete = np.arange(0, np.shape(vectors)[2])
        # for 
        # for i in range(np.shape(vectors)[1]-1, -1, -1):
        #     if i != colsToKeep[-1]:
        #         np.delete(vectors, i, 1)
        #         np.delete(vecPos, i, 1)
        #     else:
        #         np.delete(colsToKeep, -1)

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
            counter = 0
            print(noFigs)
            print(figRows)
            print(figCols)
            for i in range(0, figRows):
                for j in range(0, figCols):
                    if counter < noFigs:
                        # obj = axs[i].pcolormesh(rhoAxis, zAxis, diagnostic, cmap = "plasma", norm=LogNorm())
                        axs[i, j].pcolormesh(rhoAxis, zAxis, diagnostic[counter, :, :], cmap = "plasma", norm=LogNorm())
                        axs[i, j].quiver(vecPos2[counter, :, :, 0], vecPos2[counter, :, :, 1], vectors2[counter, :, :, 0], vectors2[counter, :, :, 1])
                        # plt.colorbar(obj)
                        axs[i, j].set_xlabel("rho")
                        axs[i, j].set_ylabel("z")
                        titleString = "Phi = " + str(phiArrayDeg[counter]) + "º"
                        axs[i, j].set_title(titleString)
                        axs[i, j].set_aspect("equal") 
                        counter += 1
        
        return
        

     

    def getLShellB(self, L, NTheta = 100, NPhi = 361):
        #Returns an array of the magnetic field at a particular L-shell 
        #deltaTheta in degrees
        #deltaPhi in degrees

        #First get the arrays set up
        #Calculate cut-off angle
        cutOffTheta = np.arcsin(np.sqrt(1/L))
        thetaArray = np.linspace(cutOffTheta, np.pi - cutOffTheta, NTheta)
        deltaPhi = 2*np.pi/(NPhi-1)

        #There are Ny rows, Nx columns, and each point has 3 things, r theta and phi
        rF_spherical = np.zeros((NPhi, NTheta, 3))
        rF_cartesian = np.zeros(np.shape(rF_spherical))
        bF = np.zeros(np.shape(rF_cartesian))
        for i in range(0, NPhi):
            for j in range(0, NTheta):
                #set theta
                rF_spherical[i, j, 1] = np.copy(thetaArray[j])
                #set phi
                rF_spherical[i, j, 2] = deltaPhi * i
                #set r
                rF_spherical[i, j, 0] = self.a*L*np.power(np.sin(rF_spherical[i, j, 1]), 2)

                 #Now have the coordinates in field-aligned spherical polars
                #Need to convert to field aligned cartesian

                theta = np.copy(rF_spherical[i, j, 1])
                phi = np.copy(rF_spherical[i, j, 2])
                rF_spherical[i, j, 1] = 0
                rF_spherical[i, j, 2] = 0
                rF_cartesian[i, j, :] = self.convertPolarToCartesian(rF_spherical[i, j, :], theta, phi)

                #Now get B field at each point
                bF[i, j, :] = self.getField(rF_cartesian[i, j, :])

                rF_spherical[i, j, 1] = theta
                rF_spherical[i, j, 2] = phi
 
               
        return rF_spherical, rF_cartesian, bF

    
    def plotDeviationColourMapLShell(self, L, NTheta = 100, NPhi = 361):
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

        thetaAxis = rF_spherical[0, :, 1]*180/np.pi
        phiAxis = rF_spherical[:, 0, 2]*180/np.pi
        # phiAxis = np.arange(0, 360 + deltaPhi, deltaPhi)
        
        ax6 = plt.figure().add_subplot()
        # imshowObject = ax6.imshow(np.transpose(ratio), cmap = "plasma", norm=LogNorm(vmin=0.1, vmax=1))
        # imshowObject = ax6.imshow(np.transpose(ratio), cmap = "plasma", norm=LogNorm())
        obj = ax6.pcolormesh(phiAxis, thetaAxis, np.transpose(ratio), cmap = "plasma", norm=LogNorm())
        ax6.set_ylim(ax6.get_ylim()[::-1])
        plt.colorbar(obj)
        ax6.set_xlabel("Phi (º)")
        ax6.set_ylabel("Theta (º)")
        titleString = "Diagnostic on L = " + str(L)
        ax6.set_title(titleString)



        return





class OTDField(Field):
    def __init__(self, m):
        self.m = m

    def getField(self, rvec):
        r = np.linalg.norm(rvec)
        return -(1/np.power(r, 3))*self.m + (3/np.power(r, 5))*np.dot(self.m, rvec)*rvec


class UniformField(Field):
    def __init__(self, fieldVec):
        self.B = fieldVec
    
    def getField(self):
        return self.B

