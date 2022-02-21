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
species = 'Proton' #Proton or Electron
arguments = np.array([q, m_p], dtype = object)
mode0 = 1 #1 for dipole, 2 for full field

L_shell = 7
phi_in = 200 * np.pi / 180
theta_in = 40 * np.pi / 180
phase = 90 * np.pi / 180
t0 = 0.
E = 1e6 * abs(q)
n = 5000000 #number of steps

#-------------------------#
if mode0 == 1:
    shape = 'Dipole'

elif mode0 == 2:
    shape = 'Quadrupole'

if species == 'Proton':
    arguments = np.array([q, m_p], dtype = object)
    
elif species == 'Electron':
    arguments = np.array([-q, m_e], dtype = object)
    
#-------------------------#

lambda_lat = (np.pi / 2) - theta_in

#Pitch angle
alpha_eq = np.arcsin(np.sqrt((np.cos(lambda_lat) ** 6) / \
                np.sqrt(1 + (3 * np.sin(lambda_lat) * np.sin(lambda_lat)))))

#Magnitude of v perpendicular
perp_mag = np.tan(alpha_eq)

#direction of zero phase
zero_phase_dir = np.array([np.cos(phi_in), np.sin(phi_in), 0])

#Rotation to determine the direction at the indicated phase
Rphase = np.array([[np.cos(phase), - np.sin(phase), 0], 
                   [np.sin(phase), np.cos(phase),   0], 
                   [0,             0,               1]])
#phase direction
phase_dir = np.matmul(Rphase, zero_phase_dir)
#initial guiding centre location
gc0_in = np.array(Sph_to_Cart(L_shell * a, np.pi / 2, phi_in))
#velocity magnitude
v_mag_in = E_to_v(E, arguments[1])
#perpendicular velocity magnitude
v0_perpendicular = np.sin(alpha_eq) * v_mag_in

B_in = B_rot_fun(gc0_in, a, g, h_coeff, mode0, R)

gyroradius0 = (gamma(v_mag_in) * arguments[1] * v0_perpendicular) / (abs(arguments[0]) * np.linalg.norm(B_in))
#Determine position of particle based on gyroradius and phase direction
r0 = gc0_in + (gyroradius0 * phase_dir)
#r0 = np.array([6., 0., 0.]) * a

#direction of velocity based on particle species
if arguments[1] == m_e:
    v_dir = np.cross(phase_dir, np.array([0, 0, 1]))
    
elif arguments[1] == m_p:
    v_dir = np.cross(np.array([0, 0, 1]), phase_dir)

direction0 = perp_mag * v_dir
direction0[2] = 1
#direction = np.array([1, 1, 1])

#%%
t, v, r, L, mew = RK(f_dvdt, t0, E, direction0, r0, n, arguments, mode0, 100)
#t, v, r, L, mew, gyroradius, gc = RK_single(f_dvdt, t0, E, direction0, r0, n, arguments, mode0, 50)

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

u_grid, v_grid = np.mgrid[0:2 * np.pi:40j, 0:np.pi:20j]
xp = np.cos(u_grid)*np.sin(v_grid)
yp = np.sin(u_grid)*np.sin(v_grid)
zp = np.cos(v_grid)

ax.plot_wireframe(xp, yp, zp, color = 'blue')

ax.plot3D(x / a, y / a, z / a, color = 'red')
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
saveRegionData('Output/RegionTests/regionTest_Uranus_7-30-200.npz', 0, '0', mode0 - 1, E / q, alpha_eq, phase, np.linalg.norm(gc0_in), np.linalg.norm(gc), gyroradius0, gyroradius)

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
gc_test = []
for i in tqdm(range(len(r))):
    current_B = B_rot_fun(r[i], a, g, h_coeff, mode0, R)
    v_par = np.dot(current_B / np.linalg.norm(current_B), v[i]) * (current_B / np.linalg.norm(current_B))
    v_perp = v[i] - v_par
    v_perp_mag = np.linalg.norm(v_perp)
    gyror = (gamma(v_perp_mag) * arguments[1] * v_perp_mag) / (abs(q) * np.linalg.norm(current_B))
    perp_vec = np.cross(v_perp, current_B)
    perp_dir = perp_vec / np.linalg.norm(perp_vec)
    gc_term = r[i] + ((arguments[0] / q) * gyror * perp_dir)
    
    gc_test.append(gc_term)
    xt.append(r[i][0])
    yt.append(r[i][1])
    zt.append(r[i][2])
    
    test.append(gyror)

gc_test = np.array(gc_test)
xt = np.array(xt)
yt = np.array(yt)
zt = np.array(zt)

gc_x=[]
gc_y=[]
gc_z=[]

for i in range(len(gc_test)):
    gc_x.append(gc_test[i][0])
    gc_y.append(gc_test[i][1])
    gc_z.append(gc_test[i][2])
    
gc_x = np.array(gc_x) 
gc_y = np.array(gc_y)
gc_z = np.array(gc_z)
#%%

plt.plot(xt[:200] / a, yt[:200] / a, color = 'blue', label = 'Particle')
plt.plot(gc_x[:200] / a, gc_y[:200] / a, 'x', color = 'red', label = 'Guiding Centre')
plt.xlabel('x ($r_U$)', fontsize=16)
plt.ylabel('y ($r_U$)', fontsize=16)
plt.title('{} {} - {:#d}keV'.format(species, shape, int(E / (q * 1e3))))
plt.legend(fontsize = 20)
plt.show()

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot3D(xt[3000:3200] / a, yt[3000:3200] / a, zt[3000:3200] / a, color = 'blue', label = 'Particle')
ax.plot3D(gc_x[3000:3200] / a, gc_y[3000:3200] / a, gc_z[3000:3200] / a, color = 'red', label = 'Guiding Centre')
ax.set_xlabel('x ($r_U$)',fontsize=16)
ax.set_ylabel('y ($r_U$)', fontsize=16)
ax.set_zlabel('z ($r_U$)', fontsize=16)
ax.set_title('{} {} - {:#d}keV'.format(species, shape, int(E / (q * 1e3))))
plt.legend(fontsize = 20)
plt.show()

#%%

savedArrays = np.load('Output/nickComparison-Proton-10000.npz')
pos = savedArrays["positions"]
vel = savedArrays["velocities"]
tim = savedArrays["times"]

#%%

data = loadRegionData('Output/RegionTests/regionTest_Uranus_7-30-200.npz')
data = selectCriteria(data, name = 'Harry', species = 'electron')

#%%
dr_dip = []
en_dip = []
dr_quad = []
en_quad = []
for i in range(len(data)):
    if data[i][3] == 0:
        dr_dip.append((data[i][8] - data[i][7]) / a)
        en_dip.append(data[i][4])
        
    if data[i][3] == 1:
        dr_quad.append((data[i][8] - data[i][7]) / a)
        en_quad.append(data[i][4])
        
# dr = (data[:, 8] - data[:, 7]) / a
# en = data[:, 4]
# gyro_plot = data[:, 10]
#plt.plot(en, dr, 'x', color = 'blue', ms = 8)
plt.plot(np.array(en_dip) / 1e3, dr_dip, 'x', color = 'red', ms = 8, label = 'Dipole')
#plt.plot(np.array(en_quad) / 1e3, dr_quad, 'x', color = 'blue', ms = 8, label = 'Quadrupole')
plt.xscale('log')
# plt.yscale('log')
titleString = "Change in equatorial r/a against initial KE for location L=" + str(np.round(L_shell)) + ", " + r"$\theta$ = " + str(np.round(theta_in * 180 / np.pi)) + ", " + r"$\phi$ = " + str(np.round(phi_in * 180 / np.pi))
plt.title(titleString, fontsize = 16)
plt.xlabel("Kinetic energy (keV)", fontsize = 16)
plt.ylabel("Change in r/a", fontsize = 16)
plt.legend(fontsize = 16)
plt.show()

#%%
plotRChangeOnEnergy(data, a, L_shell, 30, 200, logEnergy = True)
plt.show()

#%%
print('v0 =', v[0])
print('v0 magnitude =', np.linalg.norm(v[0]), 'm/s')
print('r0 =', r[0] / a, 'in r_u' )
print('gyroradius =', gyroradius0, 'm')
print('phase direction =', phase_dir)
print('B_in =', np.linalg.norm(B0), 'T')
print('v0 perp =', v0_perp,'m/s')
print('pitch angle =', alpha_eq, 'rad')

#%%
print('r0 =', r[0] / a, 'in r_u' )
print('r[-1] =', r[-1] / a)
print('gc = ', gc / a)

#%%

print(r[-1] - np.array([-107270755.39616576, -63324023.61693426, -68776417.01780729]))
# print(gc)
# print(np.linalg.norm(gc))
# print(r[0])
print(v[0])

#%%
for_save = np.array([t, v, r, L, mew, gc0_in, gc], dtype = object)
np.savez('Poster_data', for_save)
#%%
print((np.linalg.norm(gc) - np.linalg.norm(gc0_in)) / a)
#print(v[0])
#vz = v[0][2]
#%%
print(r[-1] / a)

#%%

pos = Cart_to_Sph(r0[0], r0[1], r0[2])

print(pos[0] / a)
print(pos[1] * 180 / np.pi)
print(pos[2] * 180 / np.pi)

gc0_sph = Cart_to_Sph(gc0_in[0], gc0_in[1], gc0_in[2])

print(gc0_sph[0] / a)
print(gc0_sph[1] * 180 / np.pi)
print(gc0_sph[2] * 180 / np.pi)
#%%

print(phi_in * 180 / np.pi)
#%%

print(gyroradius0 / a)

print(direction0 / np.linalg.norm(direction0))

#%%
abs_r = []

for i in r:
    
    abs_r.append(np.linalg.norm(i) / a)
    
plt.plot(t, abs_r)
plt.xlabel("time (s)", fontsize = 16)
plt.ylabel("Absolute r", fontsize = 16)

#%%
print(mew[:10])
print(0.5 * gamma(v[1]) * v[1] * v[1] / np.linalg.norm(B_rot_fun(r[1], a, g, h_coeff, 1, R)))