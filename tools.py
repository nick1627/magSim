"""
Contains code that saves and loads data so we can compare results etc.
"""
import numpy as np
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
