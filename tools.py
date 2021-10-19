"""
Contains code that saves and loads data so we can compare results etc.
"""
import numpy as np

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
