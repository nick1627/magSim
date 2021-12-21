# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 01:07:19 2021

@author: USER
"""
import numpy as np
import matplotlib.pyplot as plt
from Harry.Functions import *

def fc_k1(func, t, v, r, args = None):
    term = func(t, v, r, args)
    return term

def fc_k2(func, t, v, r, h, k1, l1, args = None):
    term = func(t + h, v + (h * k1), r + (h * l1), args)
    return term

def fc_k3(func, t, v, r, h, k1, k2, l1, l2, args = None):
    term = func(t + (h/2), v + (h * (((3 * k1) + k2) / 8)), r + (h * (((3 * l1) + l2) / 8)), args)
    return term
    
def fc_k4(func, t, v, r, h, k1, k2, k3, l1, l2, l3, args = None):
    term = func(t + (2*h/3), v + (h * (((8 * k1) + (2 * k2) + (8 * k3)) / 27)),\
                r + (h * (((8 * l1) + (2 * l2) + (8 * l3)) / 27)), args)
    return term

def fc_k5(func, t, v, r, h, k1, k2, k3, k4, l1, l2, l3, l4, args = None):
    t_term = t + (((7 - np.sqrt(21)) * h) / 14)
    v_term = v + (h * (((((9 * np.sqrt(21)) - 21) * k1) - ((8 * (7 - np.sqrt(21))) * k2)\
        + ((48 * (7 - np.sqrt(21))) * k3) - ((3 * (21 - np.sqrt(21))) * k4)) / 392))
    r_term = r + (h * (((((9 * np.sqrt(21)) - 21) * l1) - ((8 * (7 - np.sqrt(21))) * l2)\
        + ((48 * (7 - np.sqrt(21))) * l3) - ((3 * (21 - np.sqrt(21))) * l4)) / 392))
    
    term = func(t_term, v_term, r_term, args)
    return term

def fc_k6(func, t, v, r, h, k1, k2, k3, k4, k5, l1, l2, l3, l4, l5, args = None):
    t_term = t + (((7 + np.sqrt(21)) * h) / 14)
    v_term = v + (h * ((((-5 * (231 + (51 * np.sqrt(21)))) * k1) - ((40 * (7 + np.sqrt(21))) * k2)\
        - (320 * np.sqrt(21) * k3) + (3 * (21 + (121 * np.sqrt(21))) * k4) \
            + (392 * (6 + np.sqrt(21)) * k5)) / 1960))
    r_term = r + (h * ((((-5 * (231 + (51 * np.sqrt(21)))) * l1) - ((40 * (7 + np.sqrt(21))) * l2)\
        - (320 * np.sqrt(21) * l3) + (3 * (21 + (121 * np.sqrt(21))) * l4) \
            + (392 * (6 + np.sqrt(21)) * l5)) / 1960))
        
    term = func(t_term, v_term, r_term, args)
    return term

def fc_k7(func, t, v, r, h, k1, k2, k3, k4, k5, k6, l1, l2, l3, l4, l5, l6, args = None):
    t_term = t + h
    v_term = v + (h * ((((15 * (22 + (7 * np.sqrt(21)))) * k1) + (120 * k2)\
        + (40 * ((7 * np.sqrt(21)) - 5) * k3) - (63 * ((3 * np.sqrt(21)) - 2) * k4) \
            - (14 * (49 + (9 * np.sqrt(21))) * k5)) + (70 * (7 - np.sqrt(21)) * k6) / 180))
    r_term = r + (h * ((((15 * (22 + (7 * np.sqrt(21)))) * l1) + (120 * l2)\
        + (40 * ((7 * np.sqrt(21)) - 5) * l3) - (63 * ((3 * np.sqrt(21)) - 2) * l4) \
            - (14 * (49 + (9 * np.sqrt(21))) * l5)) + (70 * (7 - np.sqrt(21)) * l6) / 180))
        
    term = func(t_term, v_term, r_term, args)
    return term

def RK_coupled(f, g, h, t0, E, direction, r0, n, args):
    
    t = np.zeros(n+1)
    t[0] = t0
    
    d_norm  = direction / np.linalg.norm(direction)
    
    v_mag = E_to_v(E, args[2])
    v0 = v_mag * d_norm
    
    v = np.zeros(n+1, dtype = object)
    v[0] = v0
    r = np.zeros(n+1, dtype = object)
    r[0] = r0
    
    for i in range(n):
        k1 = fc_k1(f, t[i], v[i], r[i], args)
        l1 = fc_k1(g, t[i], v[i], r[i])
        
        k2 = fc_k2(f, t[i], v[i], r[i], h, k1, l1, args)
        l2 = fc_k2(g, t[i], v[i], r[i], h, k1, l1)
        
        k3 = fc_k3(f, t[i], v[i], r[i], h, k1, k2, l1, l2, args)
        l3 = fc_k3(g, t[i], v[i], r[i], h, k1, k2, l1, l2)
        
        k4 = fc_k4(f, t[i], v[i], r[i], h, k1, k2, k3, l1, l2, l3, args)
        l4 = fc_k4(g, t[i], v[i], r[i], h, k1, k2, k3, l1, l2, l3)
        
        k5 = fc_k5(f, t[i], v[i], r[i], h, k1, k2, k3, k4, l1, l2, l3, l4, args)
        l5 = fc_k5(g, t[i], v[i], r[i], h, k1, k2, k3, k4, l1, l2, l3, l4)
        
        k6 = fc_k6(f, t[i], v[i], r[i], h, k1, k2, k3, k4, k5, l1, l2, l3, l4, l5, args)
        l6 = fc_k6(g, t[i], v[i], r[i], h, k1, k2, k3, k4, k5, l1, l2, l3, l4, l5)
        
        # k7 = fc_k7(f, t[i], v[i], r[i], h, k1, k2, k3, k4, k5, k6, l1, l2, l3, l4, l5, l6, args)
        # l7 = fc_k7(g, t[i], v[i], r[i], h, k1, k2, k3, k4, k5, k6, l1, l2, l3, l4, l5, l6)
        
        t[i+1] = t[i] + h
        
        #print(k1, k2, k3, k4, k5, k6)
        
        v[i+1] = v[i] + ((h / 5) * ((16 * k1 / 27) + (6656 * k3 / 2565) + \
                          (28561 * k4 / 11286) - (9 * k5 / 10) + (2 * k6 / 11)))
        
        r[i+1] = r[i] + ((h / 5) * ((16 * l1 / 27) + (6656 * l3 / 2565) + \
                          (28561 * l4 / 11286) - (9 * l5 / 10) + (2 * l6 / 11)))
        
        # v[i+1] = v[i] + ((h / 180) * ((9 * k1) + (64 * k3) + (49 * k5) - (49 * k6) + (9 * k7)))
        
        # r[i+1] = r[i] + ((h / 180) * ((9 * l1) + (64 * l3) + (49 * l5) - (49 * l6) + (9 * l7)))
    
    return t, v, r