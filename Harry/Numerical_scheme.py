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
import sys, os
sys.path.insert(0, os.getcwd())
from tools import *

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

def f_dvdt(t, v, r, args, mode):
    #args[0] = q, and args[1] = m.
    dvdt = (args[0] * np.cross(v, B_rot_fun(r, a, g, h_coeff, mode, R))) / (gamma(v) * args[1])
    return dvdt

def E_to_v(E, m):
    v = c * np.sqrt(1 - (((m * c * c) / ((m * c * c) + E)) ** 2))
    return v

def gyroperiod(v, m, q, B):
   term = (gamma(v) * m * 2 * np.pi) / (abs(q) * np.linalg.norm(B))
   return term

def f_k1(func, t, v, r, h, args = None, mode = 1):
    term = h * func(t, v, r, args, mode)
    return term

def f_k2(func, t, v, r, h, k1, args = None, mode = 1):
    term = h * func(t + h, v + k1, r, args, mode)
    return term

def f_k3(func, t, v, r, h, k1, k2, args = None, mode = 1):
    term = h * func(t + (h/2), v + (((3 * k1) + k2) / 8), r, args, mode)
    return term
    
def f_k4(func, t, v, r, h, k1, k2, k3, args = None, mode = 1):
    term = h * func(t + (2*h/3), v + (((8 * k1) + (2 * k2 ) + (8 * k3)) / 27),\
                r, args, mode)
    return term

def f_k5(func, t, v, r, h, k1, k2, k3, k4, args = None, mode = 1):
    t_term = t + (((7 - np.sqrt(21)) * h) / 14)
    v_term = v + (((3 * ((3 * np.sqrt(21)) - 7) * k1) - (8 * (7 - np.sqrt(21)) * k2) \
                  + (48 * (7 - np.sqrt(21)) * k3) - (3 * (21 - np.sqrt(21)) * k4)) / 392)
    
    term = h * func(t_term, v_term, r, args, mode)
    return term

def f_k6(func, t, v, r, h, k1, k2, k3, k4, k5, args = None, mode = 1):
    t_term = t + (((7 + np.sqrt(21)) * h) / 14)
    v_term = v + (1 * ((((-5 * (231 + (51 * np.sqrt(21)))) * k1) - ((40 * (7 + np.sqrt(21))) * k2)\
        - (320 * np.sqrt(21) * k3) + (3 * (21 + (121 * np.sqrt(21))) * k4) \
            + (392 * (6 + np.sqrt(21)) * k5)) / 1960))
        
    term = h * func(t_term, v_term, r, args, mode)
    return term

def f_k7(func, t, v, r, h, k1, k2, k3, k4, k5, k6, args = None, mode = 1):
    t_term = t + h
    v_term = v + (((15 * (22 + (7 * np.sqrt(21))) * k1) + (120 * k2) + \
                   (40 * ((7 * np.sqrt(21)) - 5) * k3) - (63 * ((3 * np.sqrt(21)) - 2) * k4) \
                       - (14 * (49 + (9 * np.sqrt(21))) * k5) + (70 * (7 - np.sqrt(21)) * k6)) / 180)
        
    term = h * func(t_term, v_term, r, args, mode)
    return term

def RK(f, t0, E, direction, r0, n, args, mode, test = None):
    
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
    
    B0 = B_rot_fun(r0, a, g, h_coeff, mode, R)
    v0_par = np.dot(B0 / np.linalg.norm(B0), v0) * (B0 / np.linalg.norm(B0))
    v0_perp = np.linalg.norm(v0 - v0_par)
    mew_term0 = 0.5 * args[1] * v0_perp * v0_perp / np.linalg.norm(B0)
    mew = [mew_term0]
    
    period = gyroperiod(v0, arguments[1], arguments[0], B0)
    
    h = period / 50
    
    for i in tqdm(range(n)):
        
        if i%50 == 0:
            temp_p = gyroperiod(v[i], args[1], args[0], B_rot_fun(r[i], a, g, h_coeff, mode, R))
            h = temp_p / 50
        
        k1 = f_k1(f, t[i], v[i], r[i], h, args, mode)
        
        k2 = f_k2(f, t[i], v[i], r[i], h, k1, args, mode)
        
        k3 = f_k3(f, t[i], v[i], r[i], h, k1, k2, args, mode)
        
        k4 = f_k4(f, t[i], v[i], r[i], h, k1, k2, k3, args, mode)
        
        k5 = f_k5(f, t[i], v[i], r[i], h, k1, k2, k3, k4, args, mode)
        
        k6 = f_k6(f, t[i], v[i], r[i], h, k1, k2, k3, k4, k5, args, mode)
        
        k7 = f_k7(f, t[i], v[i], r[i], h, k1, k2, k3, k4, k5, k6, args, mode)
        
        t[i+1] = t[i] + h
 
        v[i+1] = v[i] + ((1 / 180) * ((9 * k1) + (64 * k3) + (49 * k5) + (49 * k6) + (9 * k7)))
        
        r[i+1] = r[i] + (v[i+1] * h)
        
        r_sph, theta, phi = Cart_to_Sph(r[i+1][0], r[i+1][1], r[i+1][2])
        L_term = r_sph / (a * np.sin(theta) * np.sin(theta))
        L.append(L_term)
        
        current_B = B_rot_fun(r[i+1], a, g, h_coeff, mode, R)
        v_par = np.dot(current_B / np.linalg.norm(current_B), v[i+1]) * (current_B / np.linalg.norm(current_B))
        v_perp_vec = v[i+1] - v_par
        v_perp = np.linalg.norm(v_perp_vec)
        
        mew_term = 0.5 * args[1] * v_perp * v_perp / np.linalg.norm(current_B)
        mew.append(mew_term)
        
        if test == 'Single' and r[i+1][2] < 0:
            gyroradius = (args[1] * v_perp) / (abs(args[0]) * np.linalg.norm(current_B))
            perp_vec = np.cross(v_perp_vec, current_B)
            perp_dir = perp_vec / np.linalg.norm(perp_vec)
            
            gc = r[i+1] + ((args[0] / q) * gyroradius * perp_dir)
            
            mew = np.array(mew)
            
            t = t[:i+2]
            v = v[:i+2]
            r = r[:i+2]
            
            return t, v, r, L, mew, gyroradius, gc
        
    mew = np.array(mew)
        
    return t, v, r, L, mew

#%%
#INITIAL CONDITIONS
arguments = np.array([q, m_p], dtype = object)

L_shell = 7
phi_in = 0 * np.pi / 180
theta_in = 30 * np.pi / 180
lambda_lat = (np.pi / 2) - theta_in

alpha_eq = np.arcsin(np.sqrt((np.cos(lambda_lat) ** 6) / \
                np.sqrt(1 + (3 * np.sin(lambda_lat) * np.sin(lambda_lat)))))

dir_x = np.tan(alpha_eq)    
direction = np.array([dir_x, 0, 1])
t0 = 0.
E = 1e5 * abs(q)
#direction = np.array([1, 1, 1])

#r0 = Sph_to_Cart(L_shell * a, np.pi / 2, phi_in)
r0 = np.array([6., 0., 0.]) * a

mode = 2

n = 50000
check = None#'Single'

#%%
#-------------------------#
if mode == 1:
    shape = 'Dipole'

elif mode == 2:
    shape = 'Quadrupole'

if arguments[1] == m_p:
    species = 'Proton'
    
elif arguments[1] == m_e:
    species = 'Electron'
    
if check == 'Single':
    d_norm  = direction / np.linalg.norm(direction)
    v_mag = E_to_v(E, arguments[1])
    v0 = v_mag * d_norm
    
    B0 = B_rot_fun(r0, a, g, h_coeff, mode, R)
    v0_par = np.dot(B0 / np.linalg.norm(B0), v0) * (B0 / np.linalg.norm(B0))
    v0_perp_vec = v0 - v0_par
    v0_perp = np.linalg.norm(v0_perp_vec)
    
    gyroradius0 = (arguments[1] * v0_perp) / (abs(arguments[0]) * np.linalg.norm(B0))
    perp_vec0 = np.cross(v0_perp_vec, B0)
    perp_dir0 = perp_vec0 / np.linalg.norm(perp_vec0)
    
    gc0 = r0 + ((arguments[0] / q) * gyroradius0 * perp_dir0)
#-------------------------#

if check == 'Single':
    t, v, r, L, mew, gyroradius, gc = RK(f_dvdt, t0, E, direction, r0, n, arguments, mode, check)
else:    
    t, v, r, L, mew = RK(f_dvdt, t0, E, direction, r0, n, arguments, mode, check)

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
    
    energy.append((arguments[1] * c * c) * (gamma(v[i]) - 1) / (1e3 * abs(arguments[0])))

x = np.array(x)
y = np.array(y)
z = np.array(z)

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot3D(x / a, y / a, z / a, color = 'blue')
ax.set_xlabel('x ($r_U$)',fontsize=16)
ax.set_ylabel('y ($r_U$)', fontsize=16)
ax.set_zlabel('z ($r_U$)', fontsize=16)
ax.set_title('{} {} - {:#d}keV'.format(species, shape, int(E / (q * 1e3))))
plt.show()

# plt.plot(x / a, y / a, color = 'blue')
# plt.xlabel('x ($r_U$)', fontsize=16)
# plt.ylabel('y ($r_U$)', fontsize=16)
# plt.show()

plt.plot(t_all, energy, color = 'blue')
plt.xlabel('Time (s)', fontsize=16)
plt.ylabel('Energy (keV)', fontsize=16)
plt.title('{} {} - {:#d}keV'.format(species, shape, int(E / (q * 1e3))))
plt.show()

plt.plot(t_all, z / a, color = 'blue')
plt.xlabel('Time (s)', fontsize=16)
plt.ylabel('z ($r_U$)', fontsize=16)
plt.title('{} {} - {:#d}keV'.format(species, shape, int(E / (q * 1e3))))
plt.show()

plt.plot(t_all, L, color = 'blue')
plt.xlabel('Time (s)', fontsize=16)
plt.ylabel('L', fontsize=16)
plt.title('{} {} - {:#d}keV'.format(species, shape, int(E / (q * 1e3))))
plt.show()

plt.plot(t, mew, color = 'blue')
plt.xlabel('Time (s)', fontsize=16)
plt.ylabel('Adiabatic invariant ($Am^2$)', fontsize=16)
plt.title('{} {} - {:#d}keV'.format(species, shape, int(E / (q * 1e3))))
plt.show()

#%%
#SAVE DATA
#np.savez('Harry/Simulation_data/e1000keV_1_1_1-6_0_0', t = t, v = v, r = r, L = L, mew = mew)
saveRegionData('Output/RegionTests/regionTest_Uranus_7-30-200', 0, 1, 1, E / q, alpha_eq, np.linalg.norm(gc0), np.linalg.norm(gc), gyroradius0, gyroradius)

#%%
savedArrays = np.load('Harry/Simulation_data/e1000keV_0.1_0.1_1-6_0_0.npz', allow_pickle = True)

t = savedArrays["t"]
v = savedArrays["v"]
r = savedArrays["r"]
L = savedArrays["L"]
mew = savedArrays["mew"]
#%%
x=[]
y=[]
z=[]
energy = []
for i in range(len(r)):
    x.append(r[i][0])
    y.append(r[i][1])
    z.append(r[i][2])
    
    energy.append((arguments[1] * c * c) * (gamma(v[i]) - 1) / (1e3 * abs(arguments[0])))
x = np.array(x) / a
y = np.array(y) / a
z = np.array(z) / a
#%%
plt.plot(t, energy)
plt.xlabel('Time (s)', fontsize=16)
plt.ylabel('Energy (keV)', fontsize=16)
plt.title('Dipole - 10keV')
plt.show()

plt.plot(t, z)
plt.xlabel('Time (s)', fontsize=16)
plt.ylabel('z ($r_U$)', fontsize=16)
plt.title('Dipole - 10keV')
plt.show()

#%%
print(t[-1] / 17)
#%%
print(((energy[-1] - energy[0]) / 55) * 100)

#%%
0.009198915800538998
#%%
test = []
xt = []
yt = []
zt = []
gc = []
for i in tqdm(range(len(r))):
    current_B = B_rot_fun(r[i], a, g, h_coeff, mode, R)
    v_par = np.dot(current_B / np.linalg.norm(current_B), v[i]) * (current_B / np.linalg.norm(current_B))
    v_perp = v[i] - v_par
    v_perp_mag = np.linalg.norm(v_perp)
    gyror = (arguments[1] * v_perp_mag) / (abs(q) * np.linalg.norm(current_B))
    perp_vec = np.cross(v_perp, current_B)
    perp_dir = perp_vec / np.linalg.norm(perp_vec)
    gc_term = r[i] + ((arguments[0] / q) * gyror * perp_dir)
    
    gc.append(gc_term)
    xt.append(r[i][0])
    yt.append(r[i][1])
    zt.append(r[i][2])
    
    test.append(gyror)

gc = np.array(gc)
xt = np.array(xt)
yt = np.array(yt)
zt = np.array(zt)

gc_x=[]
gc_y=[]
gc_z=[]

for i in range(len(gc)):
    gc_x.append(gc[i][0])
    gc_y.append(gc[i][1])
    gc_z.append(gc[i][2])
    
gc_x = np.array(gc_x) 
gc_y = np.array(gc_y)
gc_z = np.array(gc_z)
#%%
# plt.plot(t, test)
# plt.xlabel('Time (s)', fontsize=16)
# plt.ylabel('Gyroradius (m)', fontsize=16)
# plt.show()

# plt.plot(t, mew)
# plt.xlabel('Time (s)', fontsize=16)
# plt.ylabel('Adiabatic invariant ($Am^2$)', fontsize=16)

plt.plot(xt[:100] / a, yt[:100] / a, color = 'blue')
plt.plot(gc_x[:100] / a, gc_y[:100] / a, 'x', color = 'red')
plt.xlabel('x ($r_U$)', fontsize=16)
plt.ylabel('y ($r_U$)', fontsize=16)
plt.title('{} {} - {:#d}keV'.format(species, shape, int(E / (q * 1e3))))
plt.show()

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot3D(xt[:100] / a, yt[:100] / a, zt[:100] / a, color = 'blue')
ax.plot3D(gc_x[:100] / a, gc_y[:100] / a, gc_z[:100] / a, color = 'red')
ax.set_xlabel('x ($r_U$)',fontsize=16)
ax.set_ylabel('y ($r_U$)', fontsize=16)
ax.set_zlabel('z ($r_U$)', fontsize=16)
ax.set_title('{} {} - {:#d}keV'.format(species, shape, int(E / (q * 1e3))))
plt.show()

#%%

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
print(E/q)
#%%
print(Cart_to_Sph(r0[0], r0[1], r0[2]))