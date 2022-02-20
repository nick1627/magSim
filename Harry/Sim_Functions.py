# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 10:12:06 2022

@author: Charalambos Ioannou
14/02/2022
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


def gamma(v):
    v_mag_gam = np.linalg.norm(v)
    vc  = v_mag_gam / c
    gam = (1 - (vc * vc)) ** (-0.5)
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

def RK(f, t0, E, direction, r0, n, args, mode, step_size):
    
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
    mew_term0 = 0.5 * gamma(v_mag) * args[1] * v0_perp * v0_perp / np.linalg.norm(B0)
    mew = [mew_term0]
    
    period = gyroperiod(v0, args[1], args[0], B0)
    
    h = period / step_size
    
    for i in tqdm(range(n)):
        
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
        
        mew_term = 0.5 * gamma(v[i+1]) * args[1] * v_perp * v_perp / np.linalg.norm(current_B)
        mew.append(mew_term)
        
        if i%step_size == 0 and i > 0:
            temp_p = gyroperiod(v[i], args[1], args[0], B_rot_fun(r[i], a, g, h_coeff, mode, R))
            h = temp_p / step_size
        
        if np.linalg.norm(r[i+1]) < a :
            print('Particle has hit the planet!')
            mew = np.array(mew)
            
            t = t[:i+2]
            v = v[:i+2]
            r = r[:i+2]
            
            return t, v, r, L, mew
        
    mew = np.array(mew)
        
    return t, v, r, L, mew

def RK_single(f, t0, E, direction, r0, n, args, mode, step_size):
    
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
    mew_term0 = 0.5 * gamma(v_mag) * args[1] * v0_perp * v0_perp / np.linalg.norm(B0)
    mew = [mew_term0]
    
    period = gyroperiod(v0, args[1], args[0], B0)
    
    h = period / step_size
    
    for i in tqdm(range(n)):
        
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
        
        mew_term = 0.5 * gamma(v[i+1]) * args[1] * v_perp * v_perp / np.linalg.norm(current_B)
        mew.append(mew_term)
        
        if i%step_size == 0 and i > 0:
            temp_p = gyroperiod(v[i], args[1], args[0], B_rot_fun(r[i], a, g, h_coeff, mode, R))
            h = temp_p / step_size
        
        if r[i+1][2] < 0:
            gyroradius = (gamma(v[i+1]) * args[1] * v_perp) / (abs(args[0]) * np.linalg.norm(current_B))
            perp_vec = np.cross(v_perp_vec, current_B)
            perp_dir = perp_vec / np.linalg.norm(perp_vec)
            
            gc = r[i+1] + ((args[0] / q) * gyroradius * perp_dir)
            
            mew = np.array(mew)
            
            t = t[:i+2]
            v = v[:i+2]
            r = r[:i+2]
            
            return t, v, r, L, mew, gyroradius, gc
        
        # if np.linalg.norm(r[i+1]) < a :
        #     print('Particle has hit the planet!')
        #     mew = np.array(mew)
        #     gyroradius2 = None
        #     gc2 = None
        #     t = t[:i+2]
        #     v = v[:i+2]
        #     r = r[:i+2]
            
        #     return t, v, r, L, mew, gyroradius2, gc2
    
    print('Particle did not reach the equator')
    mew = np.array(mew)
    gyroradius2 = None
    gc2 = None
    
    return t, v, r, L, mew, gyroradius2, gc2

