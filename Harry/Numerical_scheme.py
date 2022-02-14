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
from Harry.Sim_Functions import *
from tqdm import tqdm
import sys, os
sys.path.insert(0, os.getcwd())
from tools import *

#%%
#constants
c = constants.c
m_e = constants.m_e
m_p = constants.m_p
q = constants.e

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
#INITIAL CONDITIONS
arguments = np.array([q, m_p], dtype = object)
mode = 1

L_shell = 7
phi_in = 200 * np.pi / 180
theta_in = 30 * np.pi / 180
lambda_lat = (np.pi / 2) - theta_in

alpha_eq = np.arcsin(np.sqrt((np.cos(lambda_lat) ** 6) / \
                 np.sqrt(1 + (3 * np.sin(lambda_lat) * np.sin(lambda_lat)))))

perp_mag = np.tan(alpha_eq)

zero_phase_dir = np.array([np.cos(phi_in), np.sin(phi_in), 0])
phase = 180 * np.pi / 180
Rphase = np.array([[np.cos(phase), - np.sin(phase), 0], 
                   [np.sin(phase), np.cos(phase),   0], 
                   [0,             0,               1]])
phase_dir = np.matmul(Rphase, zero_phase_dir)

gc0_in = np.array(Sph_to_Cart(L_shell * a, np.pi / 2, phi_in))

t0 = 0.
E = 1e7 * abs(q)
v_mag_in = E_to_v(E, arguments[1])

v0_perp = np.sin(alpha_eq) * v_mag_in

B0 = B_rot_fun(gc0_in, a, g, h_coeff, mode, R)

gyroradius0 = (gamma(v0_perp) * arguments[1] * v0_perp) / (abs(arguments[0]) * np.linalg.norm(B0))

r0 = gc0_in + (gyroradius0 * phase_dir)
#r0 = np.array([6., 0., 0.]) * a

if arguments[1] == m_e:
    v_dir = np.cross(phase_dir, np.array([0, 0, 1]))
    
elif arguments[1] == m_p:
    v_dir = np.cross(np.array([0, 0, 1]), phase_dir)

direction = perp_mag * v_dir
direction[2] = 1
#direction = np.array([1, 1, 1])

n = 10000000

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
    
#-------------------------#

  
#t, v, r, L, mew = RK(f_dvdt, t0, E, direction, r0, n, arguments, mode, 50)
t, v, r, L, mew, gyroradius, gc = RK_single(f_dvdt, t0, E, direction, r0, n, arguments, mode, 50)

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

plt.plot(x / a, y / a, color = 'blue')
plt.xlabel('x ($r_U$)', fontsize=16)
plt.ylabel('y ($r_U$)', fontsize=16)
plt.show()

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
saveRegionData('Output/RegionTests/regionTest_Uranus_7-30-200.npz', 0, '1', mode - 1, E / q, alpha_eq, phase, np.linalg.norm(gc0), np.linalg.norm(gc), gyroradius0, gyroradius)

#%%
savedArrays = np.load('Harry/Simulation_data/e1000keV_0.1_0.1_1-6_0_0.npz', allow_pickle = True)

t = savedArrays["t"]
v = savedArrays["v"]
r = savedArrays["r"]
L = savedArrays["L"]
mew = savedArrays["mew"]

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

plt.plot(xt[:200] / a, yt[:200] / a, color = 'blue')
plt.plot(gc_x[:200] / a, gc_y[:200] / a, 'x', color = 'red')
plt.xlabel('x ($r_U$)', fontsize=16)
plt.ylabel('y ($r_U$)', fontsize=16)
plt.title('{} {} - {:#d}keV'.format(species, shape, int(E / (q * 1e3))))
plt.show()

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot3D(xt / a, yt / a, zt / a, color = 'blue')
ax.plot3D(gc_x / a, gc_y / a, gc_z / a, color = 'red')
ax.set_xlabel('x ($r_U$)',fontsize=16)
ax.set_ylabel('y ($r_U$)', fontsize=16)
ax.set_zlabel('z ($r_U$)', fontsize=16)
ax.set_title('{} {} - {:#d}keV'.format(species, shape, int(E / (q * 1e3))))
plt.show()

#%%

savedArrays = np.load('Output/nickComparison-Proton-10000.npz')
pos = savedArrays["positions"]
vel = savedArrays["velocities"]
tim = savedArrays["times"]

#%%

data = loadRegionData('Output/RegionTests/regionTest_Uranus_7-30-200.npz')
data = selectCriteria(data, species = 'proton', field = 'dipoleOnly')

plotRChangeOnEnergy(data, a, L_shell, 30, 200)
plt.xscale('log')
plt.show()

#%%
print('v0 =', v[0])
print('v0 magnitude =', np.linalg.norm(v[0]), 'm/s')
print('r0 =', r[0] / a, 'in r_u' )
print('gyroradius =', gyroradius0, 'm')
print('phase direction =', phase_dir)
print('B0 =', np.linalg.norm(B0), 'T')
print('v0 perp =', v0_perp,'m/s')
print('pitch angle =', alpha_eq, 'rad')