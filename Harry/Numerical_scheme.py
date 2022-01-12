# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 17:50:45 2021

@author: Charalambos Ioannou
10/12/2021
"""
import numpy as np
import matplotlib.pyplot as plt
from Harry.Functions import *

#%%
#constants
c = float(299792458)
m_e = float(9.10938356e-31)
m_p = float(1.6726219e-28)
q = float(1.602176634e-19)

E_per = []
E_value = []

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

#%%
def gamma(v):
    v_mag = np.linalg.norm(v)
    gam = (1 - ((v_mag * v_mag) / (c * c))) ** (-0.5)
    return gam

def f_dvdt(t, v, r, args):
    #args[0] = q, args[1] = B, and args[2] = m.
    dvdt = (args[0] * np.cross(v, args[1])) / (gamma(v) * args[2])
    return dvdt

def f_dvdt_n(t, v, r, args):
    dvdt =  np.cross(v * c, args[0]) * np.sqrt(1 - (v * v))
    return dvdt

def E_to_v(E, m):
    v = c * np.sqrt(1 - (((m * c * c) / ((m * c * c) + E)) ** 2))
    return v

def f_k1(func, t, v, r, args = None):
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

def RK(f, h, t0, v0, r0, n, args):
    
    t = np.zeros(n+1)
    t[0] = t0
    v = np.zeros(n+1, dtype = object)
    v[0] = v0
    r = np.zeros(n+1, dtype = object)
    r[0] = r0
    
    for i in range(n):
        k1 = f_k1(f, t[i], v[i], r[i], args)
        
        k2 = f_k2(f, t[i], v[i], r[i], h, k1, args)
        
        k3 = f_k3(f, t[i], v[i], r[i], h, k1, k2, args)
        
        k4 = f_k4(f, t[i], v[i], r[i], h, k1, k2, k3, args)
        
        k5 = f_k5(f, t[i], v[i], r[i], h, k1, k2, k3, k4, args)
        
        k6 = f_k6(f, t[i], v[i], r[i], h, k1, k2, k3, k4, k5, args)
        
        k7 = f_k7(f, t[i], v[i], r[i], h, k1, k2, k3, k4, k5, k6, args)
        
        t[i+1] = t[i] + h
        
        v[i+1] = v[i] + ((1 / 180) * ((9 * k1) + (64 * k3) + (49 * k5) + (49 * k6) + (9 * k7)))
        
        r[i+1] = r[i] + (v[i+1] * h)
    
    return t, v, r

#%%
q = q
B = np.array([0., 0., 1])
arguments = np.array([q, B, m_p], dtype = object)

gyroperiod = (arguments[2] * 2 * np.pi) / (abs(q) * np.linalg.norm(B))

h = gyroperiod / 50
t0 = 0.
v0 = np.array([0., 0.1 * c, 0 * c])
E = 1e7 * q
direction = np.array([0, -1, 1])

E_value.append(E / q)

#gyroradius = (arguments[2] * np.linalg.norm(v0[:2])) / (abs(q) * np.linalg.norm(B))
r0 = np.array([0., 0., 0.])

n = 1000

#arguments2 = np.array([B, gyroperiod], dtype = object)

#t, v ,r = RK_coupled(f_dvdt, f_drdt, h, t0, E, direction, r0, n, arguments)
t, v ,r = RK(f_dvdt, h, t0, v0, r0, n, arguments)

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
    
    energy.append((arguments[2] * c * c) * (gamma(v[i]) - 1))
    
fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot3D(x, y, z)
plt.show()
plt.plot(x, y)
plt.show()

plt.plot(t_all, energy)
plt.xlabel('Time')
plt.ylabel('Energy')
plt.show()

E_per.append((energy[0] - energy[-1]) / energy[0])

#%%

plt.plot(np.array(E_value) / 1e6, np.array(E_per) * 100, 'x', ms = 10, mew = 2)
plt.xlabel('E (100 keV)', fontsize = 16)
plt.ylabel('% Energy at end of simulation after 1000 steps', fontsize = 16)
plt.show()