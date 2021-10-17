"""
Contains code that saves and loads data so we can compare results etc.
"""
import numpy as np

def saveBField(x, y, z, u, v, w, filePath):
    # This saves the data for plotting a field in the provided path
    
    np.savez(filePath, x = x, y= y, z = z, u = u, v = v, w = w)
    
    return

def loadBField(filePath):
    x, y, z, u, v, w = np.load(filePath)
    return x, y, z, u, v, w