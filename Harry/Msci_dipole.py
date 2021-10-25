# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 19:21:41 2021

@author: Charalambos Ioannou
"""
import numpy as np
import matplotlib.pyplot as plt
from Harry.Functions import *
import sys, os
sys.path.insert(0, os.getcwd())
from tools import *

#%%

#for dipole n = 1
    
a = 25600000 # Uranus' radius

#Uranus' coefficients
g01 = 11278 #z
g11 = 10928 #x
h11 = -16049 #y
g02 = 0#-9648
g12 = 0#-12284
h12 = 0#6405
g22 = 0#1453
h22 = 0#4220

args = np.array([g01, g11, h11, g02, g12, h12, g22, h22])

r = np.linspace(2 * a, 3 * a, 5)
theta = np.linspace(0, np.pi , 10)
phi = np.linspace(0, 2 * np.pi, 10)

x, y, z, u, v, w = Get_B_sph(r, theta, phi, a, args, 2)

ThreeD_plot(x, y, z, u, v, w)

#%%
a = 25600000 #Uranus radius

z1 = [0]#np.linspace(-1.5 * a, 1.5 * a, 4)
y1 = np.linspace(-1.5 * a, 1.5 * a, 20)
x1 = np.linspace(-1.5 * a, 1.5 * a, 20)

g01 = 11278 #z
g11 = 10928 #x
h11 = -16049 #y
g02 = -9648
g12 = -12284
h12 = 6405
g22 = 1453
h22 = 4220

args = np.array([g01, g11, h11, g02, g12, h12, g22, h22])

x, y, z, u, v, w = Get_B_cart(x1, y1, z1, a, args, 2, True)

TwoD_plot(x, y, u, v, 'z')
            
#%%
x_mat, y_mat, z_mat, u_mat, v_mat, w_mat = reshape(x, y, z, u, v, w, 4, 7, 5)

print(x_mat)
print(y_mat)
print(z_mat)
#%%

a = 25600000 #Uranus radius

x1 = np.linspace(-1.5 * a, 1.5 * a, 4)
y1 = np.linspace(-1.5 * a, 1.5 * a, 4)
z1 = np.linspace(-1.5 * a, 1.5 * a, 4)
    
    
g01 = 11278 #z
g11 = 10928 #x
h11 = -16049 #y
g02 = 0#-9648
g12 = 0#-12284
h12 = 0#6405
g22 = 0#1453
h22 = 0#4220


g = np.array([[g01, g11], 
              [g02, g12, g22]], dtype = object)
h = np.array([[0, h11], 
              [0, h12, h22]], dtype = object)

x, y, z, u, v, w = getB_fun(x1, y1, z1, a, g, h, 2)

TwoD_plot(y, z, v, w, 'x')

#%%

a = 25600000 #Uranus radius

z1 = [0]#np.linspace(-1.5 * a, 1.5 * a, 20)
y1 = np.linspace(-3 * a, 3 * a, 20)
x1 = np.linspace(-3 * a, 3 * a, 20)

g01 = 22454.1 #z
g11 = 0#10928 #x
h11 = 0#-16049 #y
g02 = 0#-9648
g12 = 0#-12284
h12 = 0#6405
g22 = 0#1453
h22 = 0#4220

args = np.array([g01, g11, h11, g02, g12, h12, g22, h22])

x, y, z, u, v, w = Get_B_cart(x1, y1, z1, a, args, 2, True)

TwoD_plot(x, y, u, v, 'z')

#%%

a = 25600000 #Uranus radius

x1 = np.linspace(-2* a, 2 * a, 5)
z1 = np.linspace(-2 * a, 2 * a, 5)
y1 = np.linspace(-2 * a, 2 * a, 5)

g01 = 11278 #z
g11 = 10928 #x
h11 = -16049 #y
g02 = 0#-9648
g12 = 0#-12284
h12 = 0#6405
g22 = 0#1453
h22 = 0#4220

args = np.array([g01, g11, h11, g02, g12, h12, g22, h22])

x, y, z, u, v, w = Get_B_cart(x1, y1, z1, a, args, 2)

ThreD_plot(x, y, z, u, v, w)

#%%
#Calculate rotation matrix and its inverse!

z_f = np.array([10928, -16049, 11278]) / np.sqrt((10928 * 10928) + \
                                    (-16049 * -16049) + (11278 * 11278))
z_r = np.array([0, 0, 1])
  
x_f = np.cross(z_f, z_r)
x_f = x_f / np.sqrt(np.dot(x_f, x_f))

y_f = np.cross(z_f, x_f)
y_f = y_f / np.sqrt(np.dot(y_f, y_f))

for i in range(len(x_f)):
    R[i][0] = x_f[i]
    R[i][1] = y_f[i]
    R[i][2] = z_f[i]

Rinv = np.linalg.inv(R)

R = np.array([x_f, y_f, z_f])

print(np.dot(z_r, z_r))

print(np.linalg.det(R))

#print(np.cross(x_f, y_f))
#print(np.matmul(R, np.identity(1)))
#%%

a = 25600000 #Uranus radius

z1 = [0]#np.linspace(-1.5 * a, 1.5 * a, 20)
x1 = np.linspace(-1.5 * a, 1.5 * a, 20)
y1 = np.linspace(-1.5 * a, 1.5 * a, 20)

g01 = 1#11278 #z
g11 = 0#10928 #x
h11 = 0#-16049 #y
g02 = 0#-9648
g12 = 0#-12284
h12 = 0#6405
g22 = 0#1453
h22 = 0#4220

args = np.array([g01, g11, h11, g02, g12, h12, g22, h22])

x, y, z, u, v, w = Get_B_cart_rot(x1, y1, z1, a, args, 2, np.identity(3), True)

TwoD_plot(x, y, u, v, 'z')

#%%

xn, yn, zn, un, vn, wn = loadBField('Output/dipole_nick_fieldaligned.npz')

xc, yc, zc, uc, vc, wc = reshapeCN(x, y, z, u, v, w, 4, 4, 4)

#%%

print(un - uc)

#print(x, y, z, u, v, w)