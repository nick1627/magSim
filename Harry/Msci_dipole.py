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
g01 = 1#11893 #z
g11 = 1#11579 #x
h11 = 1#- 15684 #y
g02 = 1
g12 = 1
h12 = 1
g22 = 1
h22 = 1

args = np.array([g01, g11, h11, g02, g12, h12, g22, h22])

r = [2*a]#np.linspace(2 * a, 3 * a, 1)
theta = [np.pi / 4]#np.linspace(0, np.pi, 25)
phi = [np.pi / 2]#np.linspace(0, 2 * np.pi, 25)

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

x1 = np.linspace(-1.5 * a, 1.5 * a, 4)
z1 = np.linspace(-1.5 * a, 1.5 * a, 5)
y1 = np.linspace(-1.5 * a, 1.5 * a, 7)

# r = 2*a#np.linspace(2 * a, 3 * a, 1)
# theta = np.pi / 4#np.linspace(0, np.pi, 25)
# phi = np.pi / 2#np.linspace(0, 2 * np.pi, 25)

# x, y, z = Sph_to_Cart(r, theta, phi)

# x = [x]
# y = [y]
# z = [z]

g01 = 1#11278 #z
g11 = 1#10928 #x
h11 = 1#-16049 #y
g02 = 1#-9648
g12 = 1#-12284
h12 = 1#6405
g22 = 1#1453
h22 = 1#4220

args = np.array([g01, g11, h11, g02, g12, h12, g22, h22])

x, y, z, u, v, w = Get_B_cart(x1, y1, z1, a, args, 2)

x = np.array(x) / a
y = np.array(y) / a
z = np.array(z) / a
u = np.array(u)
v = np.array(v)
w = np.array(w)

Circle1 = plt.Circle((0, 0), 1, color = 'red', fill = False)
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

# fig, ax = plt.subplots()
# ax.add_patch(Circle1)
# ax.quiver(x, y, u, v, pivot = 'mid')
# plt.title('at $z/a$ = {}'.format(0))
# plt.xlabel('x/a')
# plt.ylabel('y/a')
# #plt.savefig('h22_z=0')
plt.show()

#%%

def reshape(x, y, z, u, v, w, x_s, y_s, z_s):

    x_mat = x.reshape((x_s,y_s,z_s))
    y_mat = y.reshape((x_s,y_s,z_s))
    z_mat = z.reshape((x_s,y_s,z_s))
    u_mat = u.reshape((x_s,y_s,z_s))
    v_mat = v.reshape((x_s,y_s,z_s))
    w_mat = w.reshape((x_s,y_s,z_s))

    return x_mat, y_mat, z_mat, u_mat, v_mat, w_mat
            
#%%
x_mat, y_mat, z_mat, u_mat, v_mat, w_mat = reshape(x, y, z, u, v, w, 4, 7, 5)

print(x_mat)
print(y_mat)
print(z_mat)
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