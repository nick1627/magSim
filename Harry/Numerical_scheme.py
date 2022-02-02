# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 17:50:45 2021

@author: Charalambos Ioannou
10/12/2021
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as constants
from Harry.Functions import *
from tqdm import tqdm

#%%
a = 25600000 # Uranus' radius

#Calculate rotation matrix and its inverse!

z_f = np.array([10928, -16049, 11278]) / np.sqrt((10928 * 10928) + \
                                    (-16049 * -16049) + (11278 * 11278))
z_r = np.array([0, 0, 1])
  
x_f = np.cross(z_f, z_r)
x_f = x_f / np.sqrt(np.dot(x_f, x_f))

y_f = np.cross(z_f, x_f)
y_f = y_f / np.sqrt(np.dot(y_f, y_f))

R = np.zeros((3,3))

for i in range(len(x_f)):
    R[i][0] = x_f[i]
    R[i][1] = y_f[i]
    R[i][2] = z_f[i]

R = np.array([x_f, y_f, z_f])

Rinv = np.linalg.inv(R)

#Field coefficients

g01 = 11278 * 1e-9 #z
g11 = 10928 * 1e-9 #x
h11 = -16049 * 1e-9 #y
g02 = -9648 * 1e-9
g12 = -12284 * 1e-9
h12 = 6405 * 1e-9
g22 = 1453 * 1e-9
h22 = 4220 * 1e-9

g = np.array([[g01, g11], 
              [g02, g12, g22]], dtype = object)
h_coeff = np.array([[0, h11], 
              [0, h12, h22]], dtype = object)

#constants
c = constants.c
m_e = constants.m_e
m_p = constants.m_p
q = constants.e

# E_per = []
# E_value = []

plt.rcParams["figure.autolayout"] = True
params = {
'axes.labelsize': 12,
'font.size': 12,
#'font.family': 'sans-serif', # Optionally change the font family to sans-serif
#'font.serif': 'Arial', # Optionally change the font to Arial
'legend.fontsize': 11,
'xtick.labelsize': 12,
'ytick.labelsize': 12, 
'figure.figsize': [8, 8]
} 
plt.rcParams.update(params)

# B_mag = B_rot_fun(np.array([6 * a, 0, 0]), a, g, h_coeff, 1, R)
# #B = np.array([0, 0, np.linalg.norm(B_mag)])
# B = B_mag
# #print(np.linalg.norm(B))
#%%
def gamma(v):
    v_mag = np.linalg.norm(v)
    gam = (1 - ((v_mag * v_mag) / (c * c))) ** (-0.5)
    return gam

def f_dvdt(t, v, r, args):
    #args[0] = q, and args[1] = m.
    dvdt = (args[0] * np.cross(v, B_rot_fun(r, a, g, h_coeff, 1, R))) / (gamma(v) * args[1])
    return dvdt

def f_dvdt_n(t, v, r, args):
    dvdt =  np.cross(v * c, args[0]) * np.sqrt(1 - (v * v))
    return dvdt

def E_to_v(E, m):
    v = c * np.sqrt(1 - (((m * c * c) / ((m * c * c) + E)) ** 2))
    return v

def gyroperiod(v, m, q, B):
   term = (gamma(v) * m * 2 * np.pi) / (abs(q) * np.linalg.norm(B))
   return term

def f_k1(func, t, v, r, h, args = None):
    term = h * func(t, v, r, args)
    return term

def f_k2(func, t, v, r, h, k1, args = None):
    term = h * func(t + h, v + k1, r, args)
    return term

def f_k3(func, t, v, r, h, k1, k2, args = None):
    term = h * func(t + (h/2), v + (((3 * k1) + k2) / 8), r, args)
    return term
    
def f_k4(func, t, v, r, h, k1, k2, k3, args = None):
    term = h * func(t + (2*h/3), v + (((8 * k1) + (2 * k2 ) + (8 * k3)) / 27),\
                r, args)
    return term

def f_k5(func, t, v, r, h, k1, k2, k3, k4, args = None):
    t_term = t + (((7 - np.sqrt(21)) * h) / 14)
    v_term = v + (((3 * ((3 * np.sqrt(21)) - 7) * k1) - (8 * (7 - np.sqrt(21)) * k2) \
                  + (48 * (7 - np.sqrt(21)) * k3) - (3 * (21 - np.sqrt(21)) * k4)) / 392)
    
    term = h * func(t_term, v_term, r, args)
    return term

def f_k6(func, t, v, r, h, k1, k2, k3, k4, k5, args = None):
    t_term = t + (((7 + np.sqrt(21)) * h) / 14)
    v_term = v + (1 * ((((-5 * (231 + (51 * np.sqrt(21)))) * k1) - ((40 * (7 + np.sqrt(21))) * k2)\
        - (320 * np.sqrt(21) * k3) + (3 * (21 + (121 * np.sqrt(21))) * k4) \
            + (392 * (6 + np.sqrt(21)) * k5)) / 1960))
        
    term = h * func(t_term, v_term, r, args)
    return term

def f_k7(func, t, v, r, h, k1, k2, k3, k4, k5, k6, args = None):
    t_term = t + h
    v_term = v + (((15 * (22 + (7 * np.sqrt(21))) * k1) + (120 * k2) + \
                   (40 * ((7 * np.sqrt(21)) - 5) * k3) - (63 * ((3 * np.sqrt(21)) - 2) * k4) \
                       - (14 * (49 + (9 * np.sqrt(21))) * k5) + (70 * (7 - np.sqrt(21)) * k6)) / 180)
        
    term = h * func(t_term, v_term, r, args)
    return term

def RK(f, per, t0, E, direction, r0, n, args):
    
    t = np.zeros(n+1)
    t[0] = t0
    v = np.zeros(n+1, dtype = object)
    d_norm  = direction / np.linalg.norm(direction)
    v_mag = E_to_v(E, args[1])
    v0 = v_mag * d_norm
    v[0] = v0
    r = np.zeros(n+1, dtype = object)
    r[0] = r0
    
    r_sph, theta, phi = Cart_to_Sph(r0[0], r0[1], r0[2])
    L0 = r_sph / (a * np.sin(theta) * np.sin(theta))
    L  = [L0]
    
    h = per / 50
    
    B0 = B_rot_fun(r0, a, g, h_coeff, 1, R)
    v0_par = np.dot(B0 / np.linalg.norm(B0), v0) * (B0 / np.linalg.norm(B0))
    #v0_par = v0_par / np.linalg.norm(v0_par)
    v0_perp = np.linalg.norm(v0 - v0_par)
    mew_term0 = 0.5 * args[1] * v0_perp * v0_perp / np.linalg.norm(B0)
    mew = [mew_term0]
    
    for i in tqdm(range(n)):
        
        if i%50 == 0:
            temp_p = gyroperiod(v[i], args[1], args[0], B_rot_fun(r[i], a, g, h_coeff, 1, R))
            h = temp_p / 50
            #print(h)
        
        k1 = f_k1(f, t[i], v[i], r[i], h, args)
        
        k2 = f_k2(f, t[i], v[i], r[i], h, k1, args)
        
        k3 = f_k3(f, t[i], v[i], r[i], h, k1, k2, args)
        
        k4 = f_k4(f, t[i], v[i], r[i], h, k1, k2, k3, args)
        
        k5 = f_k5(f, t[i], v[i], r[i], h, k1, k2, k3, k4, args)
        
        k6 = f_k6(f, t[i], v[i], r[i], h, k1, k2, k3, k4, k5, args)
        
        k7 = f_k7(f, t[i], v[i], r[i], h, k1, k2, k3, k4, k5, k6, args)
        
        t[i+1] = t[i] + h
 
        v[i+1] = v[i] + ((1 / 180) * ((9 * k1) + (64 * k3) + (49 * k5) + (49 * k6) + (9 * k7)))
        
        r[i+1] = r[i] + (v[i+1] * h)
        
        r_sph, theta, phi = Cart_to_Sph(r[i+1][0], r[i+1][1], r[i+1][2])
        L_term = r_sph / (a * np.sin(theta) * np.sin(theta))
        L.append(L_term)
        
        current_B = B_rot_fun(r[i+1], a, g, h_coeff, 1, R)
        v_par = np.dot(current_B / np.linalg.norm(current_B), v[i+1]) * (current_B / np.linalg.norm(current_B))
        v_perp = np.linalg.norm(v[i+1] - v_par)
        
        mew_term = 0.5 * args[1] * v_perp * v_perp / np.linalg.norm(current_B)
        mew.append(mew_term)
        
    return t, v, r, L, mew

#%%
q = q
#B = np.array([0., 0., 1])
arguments = np.array([q, m_p], dtype = object)

t0 = 0.
#v0 = np.array([0., 0.1 * c, 0 * c])
E = 1e4 * q
direction = np.array([0.1, 0.1, 1])
d_norm  = direction / np.linalg.norm(direction)
v_mag = E_to_v(E, m_p)
v0 = v_mag * d_norm

#E_value.append(E / q)
r0 = np.array([6 * a, 0, 0.])
B_mag = B_rot_fun(r0, a, g, h_coeff, 1, R)
#B = np.array([0, 0, np.linalg.norm(B_mag)])
B = B_mag

period = gyroperiod(v0, arguments[1], arguments[0], B)

#h = period / 25

#gyroradius = (arguments[1] * np.linalg.norm(v0[:2])) / (abs(q) * np.linalg.norm(B))

n = 1000000

#arguments2 = np.array([B, gyroperiod], dtype = object)

#t, v ,r = RK_coupled(f_dvdt, f_drdt, h, t0, E, direction, r0, n, arguments)
t, v, r, L, mew = RK(f_dvdt, period, t0, E, direction, r0, n, arguments)

t_all = []
x = []
y = []
z = []
energy = []

for i in range(len(r)):
    x.append(r[i][0])
    y.append(r[i][1])
    z.append(r[i][2])
    t_all.append(t[i])
    
    energy.append((arguments[1] * c * c) * (gamma(v[i]) - 1) / (1e3 * q))

x = np.array(x)
y = np.array(y)
z = np.array(z)

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot3D(x / a, y / a, z / a)
ax.set_xlabel('x ($r_U$)',fontsize=16)
ax.set_ylabel('y ($r_U$)', fontsize=16)
ax.set_zlabel('z ($r_U$)', fontsize=16)
plt.show()

plt.plot(x / a, y / a)
plt.xlabel('x ($r_U$)', fontsize=16)
plt.ylabel('y ($r_U$)', fontsize=16)
plt.show()

plt.plot(t_all, energy)
plt.xlabel('Time (s)', fontsize=16)
plt.ylabel('Energy (keV)', fontsize=16)
plt.show()

plt.plot(t_all, z / a)
plt.xlabel('Time (s)', fontsize=16)
plt.ylabel('z ($r_U$)', fontsize=16)
plt.show()

plt.plot(t_all, L)
plt.xlabel('Time (s)', fontsize=16)
plt.ylabel('L', fontsize=16)
plt.show()
# E_per.append((energy[0] - energy[-1]) / energy[0])
#%%
plt.plot(t, mew)
plt.xlabel('Time (s)', fontsize=16)
plt.ylabel('Adiabatic invariant ($Am^2$)', fontsize=16)
plt.show()

#%%

# plt.plot(np.array(E_value) / 1e6, np.array(E_per) * 100, 'x', ms = 10, mew = 2)
# plt.xlabel('E (100 keV)', fontsize = 16)
# plt.ylabel('% Energy at end of simulation after 1000 steps', fontsize = 16)
# plt.show()

#%%

savedArrays = np.load('Output/nickComparison-Proton-10000.npz')
pos = savedArrays["positions"]
vel = savedArrays["velocities"]
tim = savedArrays["times"]

#%%
for i in range(20):
    print(r[i] - pos[i])
#%%
print(v[0], vel[0])
#%%
print(gamma(v0))

#%%
phi0 = np.arctan(r0[1] / r0[0])
phifin = np.arctan(r[-1][1] / r[-1][0])

print(phi0 * 180 / np.pi, phifin * 180 / np.pi)
#%%
print(np.linalg.norm(v0))

#%%
B0 = B_rot_fun(r0, a, g, h_coeff, 1, R)
v0_par = np.dot(B0 / np.linalg.norm(B0), v0) * (B0 / np.linalg.norm(B0))
v0_perp = np.linalg.norm(v0 - v0_par)
print(v0_perp)

v0_perp = np.linalg.norm(v0 - np.dot(B0 / np.linalg.norm(B0), v0) * (B0 / np.linalg.norm(B0)))
print(v0_perp)
#%%
print(r[:3])