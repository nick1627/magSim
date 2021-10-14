# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 19:21:41 2021

@author: Charalambos Ioannou
"""
import numpy as np
import matplotlib.pyplot as plt
from Functions import *

#%%

#for dipole n = 1
    
a = 25600000 # Uranus' radius

#Uranus' coefficients
g01 = 11893 #z
g11 = 11579 #x
h11 = - 15684 #y
g02 = 0
g12 = 0
h12 = 0
g22 = 0
h22 = 0

args = np.array([g01, g11, h11, g02, g12, h12, g22, h22])

r = np.linspace(a, 2 * a, 1)
theta = np.linspace(0.1, np.pi, 10)
phi = np.linspace(0, 2 * np.pi, 10)

x, y, z, u, v, w = Get_B_sph(r, theta, phi, a, args, 2)

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.gca(projection='3d')
#ax.set_xlim3d(-2, 2)
#ax.set_ylim3d(-2, 2)
#ax.set_zlim3d(-2, 2)

ax.quiver(np.array(x) / a, np.array(y) / a, np.array(z) / a, u,\
          v, w, length=0.4, normalize=True)
plt.show()

#%%
a = 25600000 #Uranus radius

x = np.array([0 * a])
y = np.linspace(-1.5 * a, 1.5 * a, 20)
z = np.linspace(-1.5 * a, 1.5 * a, 20)

g01 = 0#11893 #z
g11 = 0#11579 #x
h11 = 0#- 15684 #y
g02 = 1
g12 = 0
h12 = 0
g22 = 0
h22 = 0

args = np.array([g01, g11, h11, g02, g12, h12, g22, h22])

x, y, z, u, v, w = Get_B_cart(x, y, z, a, args, 2)

Circle1 = plt.Circle((0, 0), 1, color = 'red', fill = False)
params = {
   'axes.labelsize': 20,
   'font.size': 20,
   #'font.family': 'sans-serif', # Optionally change the font family to sans-serif
   #'font.serif': 'Arial', # Optionally change the font to Arial
   'legend.fontsize': 10,
   'xtick.labelsize': 20,
   'ytick.labelsize': 20, 
   'figure.figsize': [10, 10] # Using the golden ratio and standard column width of a journal
} 
plt.rcParams.update(params)

fig, ax = plt.subplots()
ax.add_patch(Circle1)
ax.quiver(np.array(y) / a, np.array(z) / a, v, w)
plt.title('at $x/a$ = {}'.format(0))
plt.xlabel('y/a')
plt.ylabel('z/a')
plt.show()