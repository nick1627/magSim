"""
Contains code that saves and loads data so we can compare results etc.
"""
from random import gauss
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import datetime as dt
from os.path import exists
import matplotlib.pyplot as plt
import copy

def saveBField(x, y, z, u, v, w, filePath):
    # This saves the data for plotting a field in the provided path
    
    np.savez(filePath, x = x, y = y, z = z, u = u, v = v, w = w)
    
    return

def loadBField(filePath):
    savedArrays = np.load(filePath)
    x = savedArrays["x"]
    y = savedArrays["y"]
    z = savedArrays["z"]
    u = savedArrays["u"]
    v = savedArrays["v"]
    w = savedArrays["w"]
    
    return x, y, z, u, v, w

def reshapeCN(x, y, z, u, v, w, x_s, y_s, z_s):
    #Changes output from 1d arrays to 3d arrays 

    x_mat = x.reshape((x_s,y_s,z_s))
    y_mat = y.reshape((x_s,y_s,z_s))
    z_mat = z.reshape((x_s,y_s,z_s))
    u_mat = u.reshape((x_s,y_s,z_s))
    v_mat = v.reshape((x_s,y_s,z_s))
    w_mat = w.reshape((x_s,y_s,z_s))

    return x_mat, y_mat, z_mat, u_mat, v_mat, w_mat

def reshapeNC(x, y, z, u, v, w):
    #Takes mag field in Nick's format
    #Returns magnetic field in Harry's format 

    #not sure if this method is correct... 

    N = np.size(x)
    xH = np.zeros(N)
    yH = np.zeros(N)
    zH = np.zeros(N)
    uH = np.zeros(N)
    vH = np.zeros(N)
    wH = np.zeros(N)
    n = 0
    for i in range(0, np.shape(x)[0]):
        for j in range(0, np.shape(x)[1]):
            for k in range(0, np.shape(x)[2]):
                xH[n] = x[i, j, k]
                yH[n] = y[i, j, k]
                zH[n] = z[i, j, k]
                uH[n] = u[i, j, k]
                vH[n] = v[i, j, k]
                wH[n] = w[i, j, k]

                n += 1

    #Order of elements may at this point be wrong
    return xH, yH, zH, uH, vH, wH




def saveRegionData(filePath, name, species, field, initialKE, pitchAngle, phase, initialRadius, finalRadius, initialGyroradius, finalGyroradius):
    """
    This function opens a file for the regional test, appends the array stored there and re-saves it.
    If the file does not already exist, a new one will be created.

    filePath:           The relative path to the file in which the data is saved
    name:               Your name, for record-keeping purposes.  It accepts your name in mulitple ways, but ultimately
                        Harry = 0, Nick = 1.
    species:            electron = 0, proton = 1
    field:              Dipole = 0, Full field = 1
    initialKE:          Float, eV
    pitchAngle:         Float, radians.  The initial pitch angle.
    phase:              Float, radians.  0 phase means in direction radially away from centre of planet.  Goes round in same direction as phi.
    initialRadius:      Float, m.  This is the radius of the centre of the gyromotion, which you must calculate.
    finalRadius:        Float, m.  This is the radius of the centre of the gyromotion, which you must calculate.
    initialGyroradius:  Float, m
    finalGyroradius:    Float, m
    """

    #Do some checks on the input name
    if name=="Harry" or name == "H" or name=="C" or name == "Charalambos" or name == 0:
        name = 0
    elif name=="Nick" or name == "N" or name == "Nicholas" or name==1:
        name = 1
    else:
        raise(Exception("Invalid name entered!"))

    if species == "electron" or species == "e" or species == "Electron" or species == "0":
        species = 0
    elif species == "proton" or species == "p" or species == "Proton" or species == "1":
        species = 1
    else:
        raise(Exception("Invalid species entered"))
    
    if field == "dipole" or field == "dipoleOnly" or field == 0:
        field = 0
    elif field == "complete" or field == "full" or field == "fullField" or field == 1:
        field = 1
    else:
        raise(Exception("Invalid field entered"))

    #construct the correct number to represent the date
    date = dt.datetime.now().month*100 + dt.datetime.now().day
    #date should be a 4 digit number, where the first two digits are the day and the second two the month

    newLine = np.array([[name, date, species, field, initialKE, pitchAngle, phase, initialRadius, finalRadius, initialGyroradius, finalGyroradius]])

    #branch depending on whether the file already exists
    if exists(filePath):
        #file exists, so we must first retrieve all data from it
        savedArray = np.load(filePath)
        oldData = savedArray["data"]
        newData = np.concatenate((oldData, newLine), axis = 0)
    else:
        #file does not exist already so we are making new one
        newData = newLine

    np.savez(filePath, data = newData)

    print("Appended line to file.")

    return

def loadRegionData(filePath):
    """
    Returns the array for the region data at the specified path
    """
    savedData = np.load(filePath)
    return savedData["data"]


def selectCriteria(data, name="", date="", species="", field="", KE="", pitchAngle="", phase = ""):
    """
    Returns an array with the desired criteria given as input

    name:       "Harry" or "Nick"
    species:    "proton" or "electron"
    field:      "fullField" or "dipoleOnly"
    phase:      Phase angle in degrees
    """
    if name == "Harry":
        deletionList = []
        for i in range(0, np.shape(data)[0]):
            if data[i, 0] != 0:
                deletionList.append(i)

        data = np.delete(data, deletionList, axis=0)

    elif name == "Nick":
        deletionList = []
        for i in range(0, np.shape(data)[0]):
            if data[i, 0] != 1:
                deletionList.append(i)

        data = np.delete(data, deletionList, axis=0)

    if species=="proton":
        deletionList = []
        for i in range(0, np.shape(data)[0]):
            if data[i, 2] != 1:
                deletionList.append(i)

        data = np.delete(data, deletionList, axis=0)
    elif species == "electron":
        deletionList = []
        for i in range(0, np.shape(data)[0]):
            if data[i, 2] != 0:
                deletionList.append(i)

        data = np.delete(data, deletionList, axis=0)

    if field=="dipoleOnly":
        deletionList = []
        for i in range(0, np.shape(data)[0]):
            if data[i, 3] != 0:
                deletionList.append(i)

        data = np.delete(data, deletionList, axis=0)
    elif field=="fullField":
        deletionList = []
        for i in range(0, np.shape(data)[0]):
            if data[i, 3] != 1:
                deletionList.append(i)

        data = np.delete(data, deletionList, axis=0)

    if phase != "":
        phase = np.pi*phase/180
        deletionList = []
        for i in range(0, np.shape(data)[0]):
            if data[i, 6] != phase:
                deletionList.append(i)

        data = np.delete(data, deletionList, axis=0)

    return data



def plotRChangeOnEnergy(regionArray, planetaryRadius, L, theta, phi, logEnergy=True):
    """
    This function accepts region data in the form produced by the loadRegionData function.
    It plots the change in r on energy.

    It is assumed that only the data that is to be plotted is input to this function.
    This is why the function accepts an array rather than the file path.
    """
    a = planetaryRadius

    #data rows are of the form:
    #name, date, species, field, initialKE, pitchAngle, phase, initialRadius, finalRadius, initialGyroradius, finalGyroradius

    #separate data by name, in case of difference between simulations
    nickData0 = []
    nickData90 = []
    nickData180 = []
    nickData270 = []
    harryData0 = []
    harryData90 = []
    harryData180 = []
    harryData270 = []

    for row in regionArray:
        if row[0] == 0: #it's harry's data
            if row[6] == 0:
                harryData0.append(copy.deepcopy(row))
            elif row[6] == 90*np.pi/180:
                harryData90.append(copy.deepcopy(row))
            elif row[6] == np.pi:
                harryData180.append(copy.deepcopy(row))
            elif row[6] == (270/180)*np.pi:
                harryData270.append(copy.deepcopy(row))
            else:
                raise(Exception("A phase is unmatched"))
        elif row[0] == 1: #it's nick's data
            if row[6] == 0:
                nickData0.append(copy.deepcopy(row))
            elif row[6] == 90*np.pi/180:
                nickData90.append(copy.deepcopy(row))
            elif row[6] == np.pi:
                nickData180.append(copy.deepcopy(row))
            elif row[6] == (270/180)*np.pi:
                nickData270.append(copy.deepcopy(row))
            else:
                raise(Exception("A phase is unmatched"))
        else:
            raise(Exception("Unmatched name!  Spooky..."))
    
    #get into array format 
    harryData0 = np.array(harryData0)
    harryData90 = np.array(harryData90)
    harryData180 = np.array(harryData180)
    harryData270 = np.array(harryData270)
    nickData0 = np.array(nickData0)
    nickData90 = np.array(nickData90)
    nickData180 = np.array(nickData180)
    nickData270 = np.array(nickData270)


    
    #now can plot

    ax = plt.figure().add_subplot()
    ax.plot(harryData0[:, 4], (harryData0[:,8] - harryData0[:,7])/a, color="red", label="H-0", linestyle = "None", marker = "x")
    ax.plot(harryData90[:, 4], (harryData90[:,8] - harryData90[:,7])/a, color="green", label="H-90", linestyle = "None", marker = "x")
    ax.plot(harryData180[:, 4], (harryData180[:,8] - harryData180[:,7])/a, color="blue", label="H-180", linestyle = "None", marker = "x")
    ax.plot(harryData270[:, 4], (harryData270[:,8] - harryData270[:,7])/a, color="black", label="H-270", linestyle = "None", marker = "x")
    ax.plot(nickData0[:, 4], (nickData0[:,8] - nickData0[:,7])/a, color="red", label="N-0", linestyle = "None", marker = "+")
    ax.plot(nickData90[:, 4], (nickData90[:,8] - nickData90[:,7])/a, color="green", label="N-90", linestyle = "None", marker = "+")
    ax.plot(nickData180[:, 4], (nickData180[:,8] - nickData180[:,7])/a, color="blue", label="N-180", linestyle = "None", marker = "+")
    ax.plot(nickData270[:, 4], (nickData270[:,8] - nickData270[:,7])/a, color="black", label="N-270", linestyle = "None", marker = "+")
    
    titleString = "Change in equatorial r/a against initial KE; L=" + str(np.round(L)) + ", " + r"$\theta$ = " + str(np.round(theta)) + ", " + r"$\phi$ = " + str(np.round(phi))
    ax.set_title(titleString)
    ax.set_xlabel("Kinetic energy (eV)")
    ax.set_ylabel("Final r/a - initial r/a")
    if logEnergy:
        ax.set_xscale('log')
    ax.legend()

    return

def plotRChangeOnEnergy2(regionArray, planetaryRadius, L, theta, phi, logEnergy=True, plotErrors=True, save=""):
    """
    This function accepts region data in the form produced by the loadRegionData function.
    It plots the change in r on energy.

    It is assumed that only the data that is to be plotted is input to this function.
    This is why the function accepts an array rather than the file path.
    """
    a = planetaryRadius

    #data rows are of the form:
    #name, date, species, field, initialKE, finalKE, pitchAngle, phase, initialRadius, finalRadius, initialGyroradius, finalGyroradius, positionError

    electronDipole = selectCriteria(regionArray, species="electron", field="dipoleOnly")
    electronFull = selectCriteria(regionArray, species="electron", field="fullField")
    protonDipole = selectCriteria(regionArray, species="proton", field="dipoleOnly")
    protonFull = selectCriteria(regionArray, species="proton", field="fullField")




    
    #now can plot
    ax = plt.figure().add_subplot()

    if plotErrors:
        ax.errorbar(protonDipole[:, 4], (protonDipole[:,9] - protonDipole[:,8])/a, yerr = protonDipole[:, 12]/a, label="Dipole, proton", color="red", marker = "+", linestyle="none")
        ax.errorbar(protonFull[:, 4], (protonFull[:,9] - protonFull[:,8])/a, yerr = protonFull[:, 12]/a, label="Full field, proton", color="purple", marker = "+", linestyle="none")
        ax.errorbar(electronDipole[:, 4], (electronDipole[:,9] - electronDipole[:,8])/a, yerr = electronDipole[:, 12]/a, label = "Dipole, electron", color="blue", marker = "x", linestyle="none")
        ax.errorbar(electronFull[:, 4], (electronFull[:,9] - electronFull[:,8])/a, yerr = electronFull[:, 12]/a, label="Full field, electron", color = "green", marker = "x", linestyle="none")
    else:
        ax.plot(protonDipole[:, 4], (protonDipole[:,9] - protonDipole[:,8])/a, label="Dipole, proton", color="red", marker = "+", linestyle="none")
        ax.plot(protonFull[:, 4], (protonFull[:,9] - protonFull[:,8])/a, label="Full field, proton", color="purple", marker = "+", linestyle="none")
        ax.plot(electronDipole[:, 4], (electronDipole[:,9] - electronDipole[:,8])/a, label = "Dipole, electron", color="blue", marker = "x", linestyle="none")
        ax.plot(electronFull[:, 4], (electronFull[:,9] - electronFull[:,8])/a, label="Full field, electron", color = "green", marker = "x", linestyle="none")
    
    
    titleString = "Change in equatorial r/a against initial KE; L=" + str(np.round(L)) + ", " + r"$\theta$ = " + str(np.round(theta)) + ", " + r"$\phi$ = " + str(np.round(phi))
    ax.set_title(titleString)
    ax.set_xlabel("Kinetic energy (eV)")
    ax.set_ylabel("Final r/a - initial r/a")
    if logEnergy:
        ax.set_xscale('log')
    ax.legend()

    if save != "":
        plt.savefig(save)

    return




def plotGyroradiusOnEnergy(regionArray):
    h = []
    n = []
    for row in regionArray:
        if row[0] == 0: #it's harry's data
            h.append(copy.deepcopy(row))
        elif row[0] == 1: #it's nick's data
            n.append(copy.deepcopy(row))
        else:
            raise(Exception("Unmatched name!  Spooky..."))
    
    h = np.array(h)
    n = np.array(n)

    
    ax = plt.figure().add_subplot()
    ax.plot(h[:, 4], h[:, 10], color="red", label="", linestyle = "None", marker = "x")
    ax.plot(n[:, 4], n[:, 10], color="red", label="", linestyle = "None", marker = "+")

    return


def deleteRegionData(path, index):
    """
    Deletes a row in the file.  Rows are identified by their
    array index.

    path:       Relative path to file
    index:      Index of line to delete
    """
    data = loadRegionData(path)
    data = np.delete(path, index, axis=0)

    np.savez(path, data = data)

    print("Deleted row " + str(index))

    return

def deleteLastRegionDataRow(path):
    """
    Deletes the last row in the data file

    path:       Relative path to file
    """

    deleteRegionData(path, -1)

    return


def deleteOlderThan(date, name, path):
    """
    Deletes data older than the given date for a given name

    date:       4 digit number.  5th September would be 0905 or just 905
    name:       String, the name of the person
    path:       The path to the data
    """

    data = loadRegionData(path)

    if name == "harry" or name == "Harry" or name == "H" or name == 0:
        name = 0
    elif name == "nick" or name == "Nick" or name == "N" or name == 1:
        name = 1
    else:
        raise(Exception("Name not recognised"))

    deletionList = []
    for i in range(0, np.shape(data)[0]):
        if data[i,0] == name:
            if data[i, 1] < date:
                deletionList.append(i)

    data = np.delete(data, deletionList, axis=0)

    np.savez(path, data = data)

    print("Data deleted.")

    return





def saveRegionData2(filePath, name, species, field, initialKE, finalKE, pitchAngle, phase, initialRadius, finalRadius, initialGyroradius, finalGyroradius, positionError, initialPhi):
    """
    This function opens a file for the regional test, appends the array stored there and re-saves it.
    If the file does not already exist, a new one will be created.

    filePath:           The relative path to the file in which the data is saved
    name:               Your name, for record-keeping purposes.  It accepts your name in mulitple ways, but ultimately
                        Harry = 0, Nick = 1.
    species:            electron = 0, proton = 1
    field:              Dipole = 0, Full field = 1
    initialKE:          Float, eV
    finalKE:            Float, eV
    pitchAngle:         Float, radians.  The initial pitch angle.
    phase:              Float, radians.  0 phase means in direction radially away from centre of planet.  Goes round in same direction as phi.
    initialRadius:      Float, m.  This is the radius of the centre of the gyromotion, which you must calculate.
    finalRadius:        Float, m.  This is the radius of the centre of the gyromotion, which you must calculate.
    initialGyroradius:  Float, m
    finalGyroradius:    Float, m
    positionError:      Float, m
    initialPhi:         The initial longitude of guiding centre in the field coord system.  Float, radians.
    """

    #Do some checks on the input name
    if name=="Harry" or name == "H" or name=="C" or name == "Charalambos" or name == 0:
        name = 0
    elif name=="Nick" or name == "N" or name == "Nicholas" or name==1:
        name = 1
    else:
        raise(Exception("Invalid name entered!"))

    if species == "electron" or species == "e" or species == "Electron" or species == "0":
        species = 0
    elif species == "proton" or species == "p" or species == "Proton" or species == "1":
        species = 1
    else:
        raise(Exception("Invalid species entered"))
    
    if field == "dipole" or field == "dipoleOnly" or field == 0:
        field = 0
    elif field == "complete" or field == "full" or field == "fullField" or field == 1:
        field = 1
    else:
        raise(Exception("Invalid field entered"))

    #construct the correct number to represent the date
    date = dt.datetime.now().month*100 + dt.datetime.now().day
    #date should be a 4 digit number, where the first two digits are the day and the second two the month

    newLine = np.array([[name, date, species, field, initialKE, finalKE, pitchAngle, phase, initialRadius, finalRadius, initialGyroradius, finalGyroradius, positionError, initialPhi]])

    #branch depending on whether the file already exists
    if exists(filePath):
        #file exists, so we must first retrieve all data from it
        savedArray = np.load(filePath)
        oldData = savedArray["data"]
        newData = np.concatenate((oldData, newLine), axis = 0)
    else:
        #file does not exist already so we are making new one
        newData = newLine

    np.savez(filePath, data = newData)

    print("Appended line to file.")

    return



def plotDeltaRHistogram(path, a, save=""):
    """
    Plots histogram for the random B field analysis.  You need a region data test file to do this.

    path:       Path to region test data file
    a:          Planetary radius
    """

    data = loadRegionData(path)
    deltaR = data[:, 9] - data[:, 8]

    deltaR = deltaR/a

    def gaussian(x, mean, stdDev, amp):
        return (amp/(stdDev*np.sqrt(2*np.pi)))*np.exp(-0.5*(x - mean)**2/(stdDev**2))

    
    
    


    ax = plt.figure().add_subplot()
    n, bins, patches = ax.hist(deltaR, bins=100)
    yData = n
    binWidth = bins[1] - bins[0]
    xData = bins + binWidth/2
    length = np.size(bins)
    xData = xData[0:length-1]

    popt, pcov = curve_fit(gaussian, xData, yData, p0 = [-0.08, 0.005, 1])

    xVals = np.linspace(bins[0], bins[-1], 1000)
    yVals = gaussian(xVals, popt[0], popt[1], popt[2])
    ax.plot(xVals, yVals)

    ax.set_ylabel("Frequency")
    ax.set_xlabel("Final r/a - initial r/a")
    ax.set_title("Effect of modifying B field parameters on net guiding centre position")

    if save!="":
        plt.savefig(save)

    print("standard deviation is %f ± %f" % (popt[1], np.sqrt(pcov[1, 1])))
    return


def plotCircumferenceGraphs(northRegionFile, southRegionFile, a, save=""):
    northData = loadRegionData(northRegionFile)
    southData = loadRegionData(southRegionFile)

    northData[:,13] *= 180/np.pi
    southData[:,13] *= 180/np.pi

    northData[:, 13] = northData[:, 13] % 360
    southData[:, 13] = southData[:, 13] % 360

    ax = plt.figure().add_subplot()
    ax.errorbar(northData[:, 13], (northData[:,9] - northData[:,8])/a, yerr = northData[:, 12]/a, label="Northbound proton", color="red", marker = "x", linestyle="none")
    ax.errorbar(southData[:, 13], (southData[:,9] - southData[:,8])/a, yerr = southData[:, 12]/a, label="Southbound proton", color = "blue", marker = "x", linestyle="none")
    ax.set_title("Change in equatorial L-shell after mirroring once")
    ax.set_xlabel("Magnetic longitude (º)")
    ax.set_ylabel("Change in equatorial L-shell (planetary radii)")

    if save != "":
        plt.savefig(save)

    ax.legend()