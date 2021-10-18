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

r = np.linspace(2 * a, 3 * a, 1)
theta = np.linspace(0, np.pi, 25)
phi = np.linspace(0, 2 * np.pi, 25)

x, y, z, u, v, w = Get_B_sph(r, theta, phi, a, args, 2)

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlim3d(-2, 2)
ax.set_ylim3d(-2, 2)
ax.set_zlim3d(-2, 2)

ax.quiver(np.array(x) / a, np.array(y) / a, np.array(z) / a, u,\
          v, w, length=0.4, normalize=True)
plt.show()

#%%
a = 25600000 #Uranus radius

x = np.linspace(-1.5 * a, 1.5 * a, 4)
y = np.linspace(-1.5 * a, 1.5 * a, 4)
z = np.linspace(-1.5 * a, 1.5 * a, 4)

g01 = 11278 #z
g11 = 10928 #x
h11 = -16049 #y
g02 = -9648
g12 = -12284
h12 = 6405
g22 = 1453
h22 = 4220

args = np.array([g01, g11, h11, g02, g12, h12, g22, h22])

x, y, z, u, v, w = Get_B_cart(x, y, z, a, args, 2)

x = np.array(x) / a
y = np.array(y) / a
z = np.array(z) / a
u = np.array(u)
v = np.array(v)
w = np.array(w)

#Circle1 = plt.Circle((0, 0), 1, color = 'red', fill = False)
params = {
   'axes.labelsize': 20,
   'font.size': 20,
   #'font.family': 'sans-serif', # Optionally change the font family to sans-serif
   #'font.serif': 'Arial', # Optionally change the font to Arial
   'legend.fontsize': 10,
   'xtick.labelsize': 20,
   'ytick.labelsize': 20, 
   'figure.figsize': [10, 10]
} 
plt.rcParams.update(params)

fig, ax = plt.subplots()
#ax.add_patch(Circle1)
#ax.quiver(y, z, v, w, pivot = 'mid')
plt.title('at $x/a$ = {}'.format(0))
plt.xlabel('y/a')
plt.ylabel('z/a')
#plt.savefig('h22_x=0')
plt.show()

#%%
x_mat = np.zeros((4,4,4))
y_mat = np.zeros((4,4,4))
z_mat = np.zeros((4,4,4))
u_mat = np.zeros((4,4,4))
v_mat = np.zeros((4,4,4))
w_mat = np.zeros((4,4,4))

for i in range(len(x)):
    s = int((i)/16)
    d = int((i)/4)%4
    f = i%4
    x_mat[s][d][f] = x[i]
    y_mat[s][d][f] = y[i]
    z_mat[s][d][f] = z[i]
    u_mat[s][d][f] = u[i]
    v_mat[s][d][f] = v[i]
    w_mat[s][d][f] = w[i]
            
#%%
print(u_mat)

#%%

a = 25600000 #Uranus radius

x = np.linspace(-1.5 * a, 1.5 * a, 4)
y = np.linspace(-1.5 * a, 1.5 * a, 4)
z = np.linspace(-1.5 * a, 1.5 * a, 4)

g01 = 11278 #z
g11 = 10928 #x
h11 = -16049 #y
g02 = -9648
g12 = -12284
h12 = 6405
g22 = 1453
h22 = 4220


g = np.array([[g01, g11], 
              [g02, g12, g22]], dtype = object)
h = np.array([[0, h11], 
              [0, h12, h22]], dtype = object)

x, y, z, u, v, w = getB_fun(x, y, z, a, g, h, 2)

x = np.array(x) / a
y = np.array(y) / a
z = np.array(z) / a
u = np.array(u)
v = np.array(v)
w = np.array(w)