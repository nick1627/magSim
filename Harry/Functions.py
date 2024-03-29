# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 11:40:41 2021

@author: Charalambos Ioannou
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

    
def Cart_to_Sph(x, y, z):
    """
    Transform Cartesian to Spherical harmonic coordinates

    Parameters
    ----------
    x : array
        x-coordinates.
    y : array
        y-coordinates.
    z : array
        z-coordinates.

    Returns
    -------
    r : array
        radius.
    theta : array
        inclinatiom.
    phi : array
        azimuth.

    """
    r = np.sqrt((x * x) + (y * y) + (z * z))
    theta = np.arctan2(np.sqrt((x * x) + (y * y)), z)
    phi = np.arctan2(y, x)
    
    if phi < 0:
        phi = phi + (2 * np.pi)
    
    return r, theta, phi

def Sph_to_Cart(r, theta, phi):
    """
    Transform spherical harmonic to Cartesian coordinates

    Parameters
    ----------
    r : array
        radius.
    theta : array
        inclinatiom.
    phi : array
        azimuth.

    Returns
    -------
    x : array
        x-coordinates.
    y : array
        y-coordinates.
    z : array
        z-coordinates
        
    """
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    
    return x, y, z

def Pnm(n, m, theta):
    """
    Calculates the Legendre Polynomials for given n, m and theta only
    up to order n = 2. 

    """
    if n == 1:
        if m == 0:
            return np.cos(theta)
        
        if m == 1:
            return np.sin(theta)
    
    if n == 2:
        if m == 0:
            term = (3 / 2) * ((np.cos(theta) * np.cos(theta)) - (1 / 3))
            return term
        
        if m == 1:
            term = (np.sqrt(3)) * np.cos(theta) * np.sin(theta)
            return term
        
        if m == 2:
            term = ((np.sqrt(3)) * 0.5) * np.sin(theta) * np.sin(theta)
            return term

def dPnm(n, m, theta):
    """
    Calculates the derivative wrt theta of the  Legendre Polynomials for given
    n, m and theta, only up to order n = 2. 

    """
    if n == 1:
        if m == 0:
            return - np.sin(theta)
        
        if m == 1:
            return np.cos(theta)
    
    if n == 2:
        if m == 0:
            term1 = (-3) * np.sin(theta) * np.cos(theta)
            return term1
        
        if m == 1:
            term2 = (np.sqrt(3)) * np.cos(2 * theta)
            return term2
        
        if m == 2:
            term3 = (np.sqrt(3)) * np.sin(theta) * np.cos(theta)
            return term3

def d2Pnm(n, m, theta):
    """
    Calculates the derivative wrt theta of the  Legendre Polynomials for given
    n, m and theta, only up to order n = 2. 

    """
    if n == 1:
        if m == 0:
            return - np.cos(theta)
        
        if m == 1:
            return - np.sin(theta)
    
    if n == 2:
        if m == 0:
            term1 = (-3) * np.cos(2 * theta)
            return term1
        
        if m == 1:
            term2 = -(2 * np.sqrt(3)) * np.sin(2 * theta)
            return term2
        
        if m == 2:
            term3 = (np.sqrt(3)) * np.cos(2 * theta)
            return term3
        
        
def Get_B_sph(r, theta, phi, a, args, n):
    """
    Calculates the B field given spherical harmonic coordinates for n = 1.
    The return values are in Cartesian coordinates.

    Parameters
    ----------
   r : array
        radius.
    theta : array
        inclinatiom.
    phi : array
        azimuth.
    a : float
        planet radius.
    g01 : float
        g01 coefficient.
    g11 : float
        g11 coefficient.
    h11 : float
        h11 coefficient.

    Returns
    -------
    x : array
        x-coordinates.
    y : array
        y-coordinates.
    z : array
        z-coordinates.
    u : array
        Bx.
    v : array
        By.
    w : array
        Bz.

    """
  
    x_all = []
    y_all = []
    z_all = []
    u_all = []
    v_all = []
    w_all = []
    
    if n == 1:
        for k in r:
            for i in theta:
                for j in phi:
                    r01 = (2 * ((a / k) ** (3))) * (args[0] * Pnm(1, 0, i))
                    r11 = (2 * ((a / k) ** (3))) * ((args[1] * np.cos(j)) + \
                                                    (args[0] * np.sin(j))) * Pnm(1, 1, i)
                    
                    B_r = r01 + r11
                    
                    theta01 = ((a / k) ** (3)) * (args[0] * dPnm(1, 0, i))
                    theta11 = (((a / k) ** (3))) * ((args[1] * np.cos(j)) + \
                                                    (args[2] * np.sin(j))) * dPnm(1, 1, i)
                    
                    B_theta = - (theta01 + theta11)
                    
                    phi01 = 0
                    phi11 = ((a / k) ** (3)) * ((args[1] * np.sin(j)) - (args[2] * np.cos(j))) \
                                                    * Pnm(1, 1, i)
                    
                    B_phi = (1 / np.sin(i)) * (phi01 + phi11)
                    
                    B = np.array([B_r, B_theta, B_phi])
            
                    x, y, z = Sph_to_Cart(k, i, j)        
            
                    u = (np.sin(i)*np.cos(j) * B[0]) + (np.cos(i)*np.cos(j) * B[1])\
                        + (-np.sin(j) * B[2])
                    v = (np.sin(i)*np.sin(j) * B[0]) + (np.cos(i)*np.sin(j) * B[1])\
                        +  (np.cos(j) * B[2])
                    w = (np.cos(i) * B[0]) + (-np.sin(i) * B[1])
            
                    x_all.append(x)
                    y_all.append(y)
                    z_all.append(z)
                    u_all.append(u)
                    v_all.append(v)
                    w_all.append(w)
    
    if n == 2:
        for k in r:
            for i in theta:
                for j in phi:
                    r01 = (2 * ((a / k) ** 3)) * (args[0] * Pnm(1, 0, i))
                    r11 = (2 * ((a / k) ** 3)) * ((args[1] * np.cos(j)) + \
                                                    (args[2] * np.sin(j))) * Pnm(1, 1, i)
                    r02 = (3 * ((a / k) ** 4)) * (args[3] * Pnm(2, 0, i))
                    r12 = (3 * ((a / k) ** 4)) * ((args[4] * np.cos(j)) + \
                                        (args[5] * np.sin(j))) * Pnm(2, 1, i)
                    r22 = (3 * ((a / k) ** 4)) * ((args[6] * np.cos(2 * j)) + \
                                        (args[7] * np.sin(2 * j))) * Pnm(2, 2, i)
                    
                    B_r = r01 + r11 + r02 + r12 + r22
                    
                    theta01 = ((a / k) ** (3)) * (args[0] * dPnm(1, 0, i))
                    theta11 = (((a / k) ** (3))) * ((args[1] * np.cos(j)) + \
                                                    (args[2] * np.sin(j))) * dPnm(1, 1, i)
                    theta02 = ((a / k) ** 4) * (args[3] * dPnm(2, 0, i))
                    theta12 = ((a / k) ** 4) * ((args[4] * np.cos(j)) + \
                                        (args[5] * np.sin(j))) * dPnm(2, 1, i)
                    theta22 = ((a / k) ** 4) * ((args[6] * np.cos(2 * j)) + \
                                        (args[7] * np.sin(2 * j))) * dPnm(2, 2, i)
                    
                    B_theta = - (theta01 + theta11 + theta02 + theta12 + theta22)
                    
                    phi01 = 0
                    phi11 = ((a / k) ** (3)) * ((args[1] * np.sin(j)) - (args[2] * np.cos(j))) \
                                                    * Pnm(1, 1, i)
                    phi02 = 0
                    phi12 = ((a / k) ** 4) * ((args[4] * np.sin(j)) - \
                                        (args[5] * np.cos(j))) * Pnm(2, 1, i)
                    phi22 = (2 * ((a / k) ** 4)) * ((args[6] * np.sin(2 * j)) - \
                                        (args[7] * np.cos(2 * j))) * Pnm(2, 2, i)
                    
                    B_phi = (1 / np.sin(i)) * (phi01 + phi11 + phi02 + phi12 + phi22)
                    
                    B = np.array([B_r, B_theta, B_phi])
                    
                    x, y, z = Sph_to_Cart(k, i, j)        
            
                    u = (np.sin(i)*np.cos(j) * B[0]) + (np.cos(i)*np.cos(j) * B[1])\
                        + (-np.sin(j) * B[2])
                    v = (np.sin(i)*np.sin(j) * B[0]) + (np.cos(i)*np.sin(j) * B[1])\
                        +  (np.cos(j) * B[2])
                    w = (np.cos(i) * B[0]) + (-np.sin(i) * B[1])
            
                    x_all.append(x)
                    y_all.append(y)
                    z_all.append(z)
                    u_all.append(u)
                    v_all.append(v)
                    w_all.append(w)
                    
    x_all = np.array(x_all) / a
    y_all = np.array(y_all) / a
    z_all = np.array(z_all) / a
    u_all = np.array(u_all)
    v_all = np.array(v_all)
    w_all = np.array(w_all) 
    
    return x_all, y_all, z_all, u_all, v_all, w_all

def Get_B_sph_rot(r_in, theta_in, phi_in, a, args, n, R):
    """
    Calculates the B field given spherical harmonic coordinates for n = 1.
    The return values are in Cartesian coordinates.

    Parameters
    ----------
   r : array
        radius.
    theta : array
        inclinatiom.
    phi : array
        azimuth.
    a : float
        planet radius.
    g01 : float
        g01 coefficient.
    g11 : float
        g11 coefficient.
    h11 : float
        h11 coefficient.

    Returns
    -------
    x : array
        x-coordinates.
    y : array
        y-coordinates.
    z : array
        z-coordinates.
    u : array
        Bx.
    v : array
        By.
    w : array
        Bz.

    """
  
    x_all = []
    y_all = []
    z_all = []
    u_all = []
    v_all = []
    w_all = []
    
    if n == 1:
        for i in r_in:
            for j in theta_in:
                for k in phi_in:
                    x, y, z = Sph_to_Cart(i, j, k)
                    xyz = np.array([x, y, z])
                    x_r, y_r, z_r = np.matmul(np.linalg.inv(R), xyz)
                    r, theta, phi = Cart_to_Sph(x_r, y_r, z_r)
                    
                    r01 = (2 * ((a / r) ** 3)) * (args[0] * Pnm(1, 0, theta))
                    r11 = (2 * ((a / r) ** 3)) * ((args[1] * np.cos(phi)) + \
                                        (args[2] * np.sin(phi))) * Pnm(1, 1, theta)
                    
                    B_r = r01 + r11
                    
                    theta01 = ((a / r) ** 3) * (args[0] * (dPnm(1, 0, theta)))
                    theta11 = ((a / r) ** 3) * ((args[1] * np.cos(phi)) + \
                                      (args[2] * np.sin(phi))) * dPnm(1, 1, theta)
                    
                    B_theta = - (theta01 + theta11)
                    
                    phi01 = 0
                    phi11 = ((a / r) ** (3)) * ((args[1] * np.sin(phi)) - (args[2] * np.cos(phi))) \
                                                    * Pnm(1, 1, theta)
                    
                    B_phi = (1 / np.sin(theta)) * (phi01 + phi11)
                    
                    B = np.array([B_r, B_theta, B_phi])        
            
                    u = (np.sin(theta)*np.cos(phi) * B[0]) + (np.cos(theta)*np.cos(phi) * B[1])\
                        + (-np.sin(phi) * B[2])
                    v = (np.sin(theta)*np.sin(phi) * B[0]) + (np.cos(theta)*np.sin(phi) * B[1])\
                        +  (np.cos(phi) * B[2])
                    w = (np.cos(theta) * B[0]) + (-np.sin(theta) * B[1])   
                    
                    B_rot = np.array([u, v, w])
                    B_f = np.matmul(R, B_rot)
                    
                    x_all.append(x)
                    y_all.append(y)
                    z_all.append(z)
                    u_all.append(B_f[0])
                    v_all.append(B_f[1])
                    w_all.append(B_f[2])
    
    if n == 2:
        for i in r_in:
            for j in theta_in:
                for k in phi_in:
                    x, y, z = Sph_to_Cart(i, j, k)
                    xyz = np.array([x, y, z])
                    x_r, y_r, z_r = np.matmul(np.linalg.inv(R), xyz)
                    r, theta, phi = Cart_to_Sph(x_r, y_r, z_r)
                    
                    r01 = (2 * ((a / r) ** 3)) * (args[0] * Pnm(1, 0, theta))
                    r11 = (2 * ((a / r) ** 3)) * ((args[1] * np.cos(phi)) + \
                                        (args[2] * np.sin(phi))) * Pnm(1, 1, theta)
                    r02 = (3 * ((a / r) ** 4)) * (args[3] * Pnm(2, 0, theta))
                    r12 = (3 * ((a / r) ** 4)) * ((args[4] * np.cos(phi)) + \
                                        (args[5] * np.sin(phi))) * Pnm(2, 1, theta)
                    r22 = (3 * ((a / r) ** 4)) * ((args[6] * np.cos(2 * phi)) + \
                                        (args[7] * np.sin(2 * phi))) * Pnm(2, 2, theta)
                    
                    B_r = r01 + r11 + r02 + r12 + r22
                    
                    theta01 = ((a / r) ** 3) * (args[0] * (dPnm(1, 0, theta)))
                    theta11 = ((a / r) ** 3) * ((args[1] * np.cos(phi)) + \
                                      (args[2] * np.sin(phi))) * dPnm(1, 1, theta)
                    theta02 = ((a / r) ** 4) * (args[3] * dPnm(2, 0, theta))
                    theta12 = ((a / r) ** 4) * ((args[4] * np.cos(phi)) + \
                                        (args[5] * np.sin(phi))) * dPnm(2, 1, theta)
                    theta22 = ((a / r) ** 4) * ((args[6] * np.cos(2 * phi)) + \
                                        (args[7] * np.sin(2 * phi))) * dPnm(2, 2, theta)
                    
                    B_theta = - (theta01 + theta11 + theta02 + theta12 + theta22)
                    
                    phi01 = 0
                    phi11 = ((a / r) ** (3)) * ((args[1] * np.sin(phi)) - (args[2] * np.cos(phi))) \
                                                    * Pnm(1, 1, theta)
                    phi02 = 0
                    phi12 = ((a / r) ** 4) * ((args[4] * np.sin(phi)) - \
                                        (args[5] * np.cos(phi))) * Pnm(2, 1, theta)
                    phi22 = (2 * ((a / r) ** 4)) * ((args[6] * np.sin(2 * phi)) - \
                                        (args[7] * np.cos(2 * phi))) * Pnm(2, 2, theta)
                    
                    B_phi = (1 / np.sin(theta)) * (phi01 + phi11 + phi02 + phi12 + phi22)
                    
                    B = np.array([B_r, B_theta, B_phi])       
                    
                    u = (np.sin(theta)*np.cos(phi) * B[0]) + (np.cos(theta)*np.cos(phi) * B[1])\
                        + (-np.sin(phi) * B[2])
                    v = (np.sin(theta)*np.sin(phi) * B[0]) + (np.cos(theta)*np.sin(phi) * B[1])\
                        +  (np.cos(phi) * B[2])
                    w = (np.cos(theta) * B[0]) + (-np.sin(theta) * B[1])
                    
                    B_rot = np.array([u, v, w])
                    B_f = np.matmul(R, B_rot)
                    
                    x_all.append(x)
                    y_all.append(y)
                    z_all.append(z)
                    u_all.append(B_f[0])
                    v_all.append(B_f[1])
                    w_all.append(B_f[2])
                    
    x_all = np.array(x_all) / a
    y_all = np.array(y_all) / a
    z_all = np.array(z_all) / a
    u_all = np.array(u_all)
    v_all = np.array(v_all)
    w_all = np.array(w_all)
                    
    return x_all, y_all, z_all, u_all, v_all, w_all

def Get_B_cart(x, y, z, a, args, n, planet = None):
    """
    Calculates the B field given Cartesian coordinates for n = 1.
    The return values are in Cartesian coordinates.
    Parameters
    ----------
    x : array
        x-coordinates.
    y : array
        y-coordinates.
    z : array
        z-coordinates.
    a : float
        planet radius.
    args : array.
        contains the coefficients used to calculate B.
    n : integer.
        The order of the series. Could be n=1 or n=2 only. 
    Returns
    -------
    x : array
        x-coordinates.
    y : array
        y-coordinates.
    z : array
        z-coordinates.
    u : array
        Bx.
    v : array
        By.
    w : array
        Bz.
    """
    
    x_all = []
    y_all = []
    z_all = []
    u_all = []
    v_all = []
    w_all = []
    
    if n == 1:
    
        for k in x:
            for i in y:
                for j in z:
                    r, theta, phi = Cart_to_Sph(k, i, j)
                    
                    r01 = (2 * ((a / r) ** 3)) * (args[0] * Pnm(1, 0, theta))
                    r11 = (2 * ((a / r) ** 3)) * ((args[1] * np.cos(phi)) + \
                                        (args[2] * np.sin(phi))) * Pnm(1, 1, theta)
                    
                    B_r = r01 + r11
                    
                    theta01 = ((a / r) ** 3) * (args[0] * (dPnm(1, 0, theta)))
                    theta11 = ((a / r) ** 3) * ((args[1] * np.cos(phi)) + \
                                      (args[2] * np.sin(phi))) * dPnm(1, 1, theta)
                    
                    B_theta = - (theta01 + theta11)
                    
                    phi01 = 0
                    phi11 = ((a / r) ** (3)) * ((args[1] * np.sin(phi)) - (args[2] * np.cos(phi))) \
                                                    * Pnm(1, 1, theta)
                    
                    B_phi = (1 / np.sin(theta)) * (phi01 + phi11)
                    
                    B = np.array([B_r, B_theta, B_phi])        
            
                    u = (np.sin(theta)*np.cos(phi) * B[0]) + (np.cos(theta)*np.cos(phi) * B[1])\
                        + (-np.sin(phi) * B[2])
                    v = (np.sin(theta)*np.sin(phi) * B[0]) + (np.cos(theta)*np.sin(phi) * B[1])\
                        +  (np.cos(phi) * B[2])
                    w = (np.cos(theta) * B[0]) + (-np.sin(theta) * B[1])   
                    
                    if planet == True:
                        if np.sqrt((k * k) + (i * i) + (j * j)) > a:
                
                            x_all.append(k)
                            y_all.append(i)
                            z_all.append(j)
                            u_all.append(u)
                            v_all.append(v)
                            w_all.append(w)
                
                    else:
                        x_all.append(k)
                        y_all.append(i)
                        z_all.append(j)
                        u_all.append(u)
                        v_all.append(v)
                        w_all.append(w)
    
    if n == 2:
        for k in x:
            for i in y:
                for j in z:
                    r, theta, phi = Cart_to_Sph(k, i, j)
                    
                    r01 = (2 * ((a / r) ** 3)) * (args[0] * Pnm(1, 0, theta))
                    r11 = (2 * ((a / r) ** 3)) * ((args[1] * np.cos(phi)) + \
                                        (args[2] * np.sin(phi))) * Pnm(1, 1, theta)
                    r02 = (3 * ((a / r) ** 4)) * (args[3] * Pnm(2, 0, theta))
                    r12 = (3 * ((a / r) ** 4)) * ((args[4] * np.cos(phi)) + \
                                        (args[5] * np.sin(phi))) * Pnm(2, 1, theta)
                    r22 = (3 * ((a / r) ** 4)) * ((args[6] * np.cos(2 * phi)) + \
                                        (args[7] * np.sin(2 * phi))) * Pnm(2, 2, theta)
                    
                    B_r = r01 + r11 + r02 + r12 + r22
                    
                    theta01 = ((a / r) ** 3) * (args[0] * (dPnm(1, 0, theta)))
                    theta11 = ((a / r) ** 3) * ((args[1] * np.cos(phi)) + \
                                      (args[2] * np.sin(phi))) * dPnm(1, 1, theta)
                    theta02 = ((a / r) ** 4) * (args[3] * dPnm(2, 0, theta))
                    theta12 = ((a / r) ** 4) * ((args[4] * np.cos(phi)) + \
                                        (args[5] * np.sin(phi))) * dPnm(2, 1, theta)
                    theta22 = ((a / r) ** 4) * ((args[6] * np.cos(2 * phi)) + \
                                        (args[7] * np.sin(2 * phi))) * dPnm(2, 2, theta)
                    
                    B_theta = - (theta01 + theta11 + theta02 + theta12 + theta22)
                    
                    phi01 = 0
                    phi11 = ((a / r) ** (3)) * ((args[1] * np.sin(phi)) - (args[2] * np.cos(phi))) \
                                                    * Pnm(1, 1, theta)
                    phi02 = 0
                    phi12 = ((a / r) ** 4) * ((args[4] * np.sin(phi)) - \
                                        (args[5] * np.cos(phi))) * Pnm(2, 1, theta)
                    phi22 = (2 * ((a / r) ** 4)) * ((args[6] * np.sin(2 * phi)) - \
                                        (args[7] * np.cos(2 * phi))) * Pnm(2, 2, theta)
                    
                    B_phi = (1 / np.sin(theta)) * (phi01 + phi11 + phi02 + phi12 + phi22)
                    
                    B = np.array([B_r, B_theta, B_phi])        
                    
                    u = (np.sin(theta)*np.cos(phi) * B[0]) + (np.cos(theta)*np.cos(phi) * B[1])\
                        + (-np.sin(phi) * B[2])
                    v = (np.sin(theta)*np.sin(phi) * B[0]) + (np.cos(theta)*np.sin(phi) * B[1])\
                        +  (np.cos(phi) * B[2])
                    w = (np.cos(theta) * B[0]) + (-np.sin(theta) * B[1])
                    
                    if planet == True:
                        if np.sqrt((k * k) + (i * i) + (j * j)) > a:
                
                            x_all.append(k)
                            y_all.append(i)
                            z_all.append(j)
                            u_all.append(u)
                            v_all.append(v)
                            w_all.append(w)
                
                    else:
                        x_all.append(k)
                        y_all.append(i)
                        z_all.append(j)
                        u_all.append(u)
                        v_all.append(v)
                        w_all.append(w)
    
    x_all = np.array(x_all) / a
    y_all = np.array(y_all) / a
    z_all = np.array(z_all) / a
    u_all = np.array(u_all)
    v_all = np.array(v_all)
    w_all = np.array(w_all)                
    
    return x_all, y_all, z_all, u_all, v_all, w_all

def Get_B_cart_rot(x, y, z, a, args, n, R, planet = None):
    """
    Calculates the B field given Cartesian coordinates for n = 1.
    The return values are in Cartesian coordinates.

    Parameters
    ----------
    x : array
        x-coordinates.
    y : array
        y-coordinates.
    z : array
        z-coordinates.
    a : float
        planet radius.
    args : array.
        contains the coefficients used to calculate B.
    n : integer.
        The order of the series. Could be n=1 or n=2 only. 

    Returns
    -------
    x : array
        x-coordinates.
    y : array
        y-coordinates.
    z : array
        z-coordinates.
    u : array
        Bx.
    v : array
        By.
    w : array
        Bz.

    """
    
    x_all = []
    y_all = []
    z_all = []
    u_all = []
    v_all = []
    w_all = []
    
    if n == 1:
    
        for k in x:
            for i in y:
                for j in z:
                    xyz = np.array([k, i, j])
                    #print(len(xyz))
                    x_r, y_r, z_r = np.matmul(np.linalg.inv(R), xyz)
                    r, theta, phi = Cart_to_Sph(x_r, y_r, z_r)
                    
                    r01 = (2 * ((a / r) ** 3)) * (args[0] * Pnm(1, 0, theta))
                    r11 = (2 * ((a / r) ** 3)) * ((args[1] * np.cos(phi)) + \
                                        (args[2] * np.sin(phi))) * Pnm(1, 1, theta)
                    
                    B_r = r01 + r11
                    
                    theta01 = ((a / r) ** 3) * (args[0] * (dPnm(1, 0, theta)))
                    theta11 = ((a / r) ** 3) * ((args[1] * np.cos(phi)) + \
                                      (args[2] * np.sin(phi))) * dPnm(1, 1, theta)
                    
                    B_theta = - (theta01 + theta11)
                    
                    phi01 = 0
                    phi11 = ((a / r) ** (3)) * ((args[1] * np.sin(phi)) - (args[2] * np.cos(phi))) \
                                                    * Pnm(1, 1, theta)
                    
                    B_phi = (1 / np.sin(theta)) * (phi01 + phi11)
                    
                    B = np.array([B_r, B_theta, B_phi])        
            
                    u = (np.sin(theta)*np.cos(phi) * B[0]) + (np.cos(theta)*np.cos(phi) * B[1])\
                        + (-np.sin(phi) * B[2])
                    v = (np.sin(theta)*np.sin(phi) * B[0]) + (np.cos(theta)*np.sin(phi) * B[1])\
                        +  (np.cos(phi) * B[2])
                    w = (np.cos(theta) * B[0]) + (-np.sin(theta) * B[1])   
                    
                    B_rot = np.array([u, v, w])
                    B_f = np.matmul(R, B_rot)
                    
                    if planet == True:
                    
                        if np.sqrt((k * k) + (i * i) + (j * j)) > a:
                        
                            x_all.append(k)
                            y_all.append(i)
                            z_all.append(j)
                            u_all.append(B_f[0])
                            v_all.append(B_f[1])
                            w_all.append(B_f[2])
                    
                    else:
                        x_all.append(k)
                        y_all.append(i)
                        z_all.append(j)
                        u_all.append(B_f[0])
                        v_all.append(B_f[1])
                        w_all.append(B_f[2])
    
    if n == 2:
        for k in x:
            for i in y:
                for j in z:
                    xyz = np.array([k, i, j])
                    x_r, y_r, z_r = np.matmul(np.linalg.inv(R), xyz)
                    r, theta, phi = Cart_to_Sph(x_r, y_r, z_r)
                    
                    r01 = (2 * ((a / r) ** 3)) * (args[0] * Pnm(1, 0, theta))
                    r11 = (2 * ((a / r) ** 3)) * ((args[1] * np.cos(phi)) + \
                                        (args[2] * np.sin(phi))) * Pnm(1, 1, theta)
                    r02 = (3 * ((a / r) ** 4)) * (args[3] * Pnm(2, 0, theta))
                    r12 = (3 * ((a / r) ** 4)) * ((args[4] * np.cos(phi)) + \
                                        (args[5] * np.sin(phi))) * Pnm(2, 1, theta)
                    r22 = (3 * ((a / r) ** 4)) * ((args[6] * np.cos(2 * phi)) + \
                                        (args[7] * np.sin(2 * phi))) * Pnm(2, 2, theta)
                    
                    B_r = r01 + r11 + r02 + r12 + r22
                    
                    theta01 = ((a / r) ** 3) * (args[0] * (dPnm(1, 0, theta)))
                    theta11 = ((a / r) ** 3) * ((args[1] * np.cos(phi)) + \
                                      (args[2] * np.sin(phi))) * dPnm(1, 1, theta)
                    theta02 = ((a / r) ** 4) * (args[3] * dPnm(2, 0, theta))
                    theta12 = ((a / r) ** 4) * ((args[4] * np.cos(phi)) + \
                                        (args[5] * np.sin(phi))) * dPnm(2, 1, theta)
                    theta22 = ((a / r) ** 4) * ((args[6] * np.cos(2 * phi)) + \
                                        (args[7] * np.sin(2 * phi))) * dPnm(2, 2, theta)
                    
                    B_theta = - (theta01 + theta11 + theta02 + theta12 + theta22)
                    
                    phi01 = 0
                    phi11 = ((a / r) ** (3)) * ((args[1] * np.sin(phi)) - (args[2] * np.cos(phi))) \
                                                    * Pnm(1, 1, theta)
                    phi02 = 0
                    phi12 = ((a / r) ** 4) * ((args[4] * np.sin(phi)) - \
                                        (args[5] * np.cos(phi))) * Pnm(2, 1, theta)
                    phi22 = (2 * ((a / r) ** 4)) * ((args[6] * np.sin(2 * phi)) - \
                                        (args[7] * np.cos(2 * phi))) * Pnm(2, 2, theta)
                    
                    B_phi = (1 / np.sin(theta)) * (phi01 + phi11 + phi02 + phi12 + phi22)
                    
                    B = np.array([B_r, B_theta, B_phi])       
                    
                    u = (np.sin(theta)*np.cos(phi) * B[0]) + (np.cos(theta)*np.cos(phi) * B[1])\
                        + (-np.sin(phi) * B[2])
                    v = (np.sin(theta)*np.sin(phi) * B[0]) + (np.cos(theta)*np.sin(phi) * B[1])\
                        +  (np.cos(phi) * B[2])
                    w = (np.cos(theta) * B[0]) + (-np.sin(theta) * B[1])
                    
                    B_rot = np.array([u, v, w])
                    B_f = np.matmul(R, B_rot)
                    
                    if planet == True:
                    
                        if np.sqrt((k * k) + (i * i) + (j * j)) > a:
                        
                            x_all.append(k)
                            y_all.append(i)
                            z_all.append(j)
                            u_all.append(B_f[0])
                            v_all.append(B_f[1])
                            w_all.append(B_f[2])
                    
                    else:
                        x_all.append(k)
                        y_all.append(i)
                        z_all.append(j)
                        u_all.append(B_f[0])
                        v_all.append(B_f[1])
                        w_all.append(B_f[2])
                    
    x_all = np.array(x_all) / a
    y_all = np.array(y_all) / a
    z_all = np.array(z_all) / a
    u_all = np.array(u_all)
    v_all = np.array(v_all)
    w_all = np.array(w_all)
                    
    return x_all, y_all, z_all, u_all, v_all, w_all

def B_sph(r, theta, phi, a, g, h, n):
    
    frac = a/r
    
    B_r = 0
    B_theta = 0
    B_phi = 0
    
    for l in range(1, n+1):
        frac2 = frac ** (l+2)
        
        for m in range(l+1):
            
            B_r += (l+1) * frac2 * ((g[l-1][m] * np.cos(m * phi)) \
                + (h[l-1][m] * np.sin(m * phi))) * Pnm(l, m, theta)
            
            B_theta -= frac2 * ((g[l-1][m] * np.cos(m * phi)) \
                    + (h[l-1][m] * np.sin(m * phi))) * dPnm(l, m, theta)
            
            B_phi += (1 / np.sin(theta)) * m * frac2 * ((g[l-1][m] * np.sin(m * phi)) \
                    - (h[l-1][m] * np.cos(m * phi))) * Pnm(l, m, theta)
    
    B = np.array([B_r, B_theta, B_phi])    
    
    return B
    
def getB_fun(x, y, z, a, g, h, n, planet = None):
    
    x_all = []
    y_all = []
    z_all = []
    u_all = []
    v_all = []
    w_all = []
    
    for k in x:
        for i in y:
            for j in z:
                r, theta, phi = Cart_to_Sph(k, i, j)
                
                frac = a/r
                
                B_r = 0
                B_theta = 0
                B_phi = 0
                
                for l in range(1, n+1):
                    frac2 = frac ** (l+2)
                    
                    for m in range(l+1):
                        
                        B_r += (l+1) * frac2 * ((g[l-1][m] * np.cos(m * phi)) \
                            + (h[l-1][m] * np.sin(m * phi))) * Pnm(l, m, theta)
                        
                        B_theta -= frac2 * ((g[l-1][m] * np.cos(m * phi)) \
                                + (h[l-1][m] * np.sin(m * phi))) * dPnm(l, m, theta)
                        
                        B_phi += (1 / np.sin(theta)) * m * frac2 * ((g[l-1][m] * np.sin(m * phi)) \
                                - (h[l-1][m] * np.cos(m * phi))) * Pnm(l, m, theta)
                
                B = np.array([B_r, B_theta, B_phi])        
                    
                u = (np.sin(theta)*np.cos(phi) * B[0]) + (np.cos(theta)*np.cos(phi) * B[1])\
                    + (-np.sin(phi) * B[2])
                v = (np.sin(theta)*np.sin(phi) * B[0]) + (np.cos(theta)*np.sin(phi) * B[1])\
                    +  (np.cos(phi) * B[2])
                w = (np.cos(theta) * B[0]) + (-np.sin(theta) * B[1])
                
                if planet == True:
                    if np.sqrt((k * k) + (i * i) + (j * j)) > a:
                
                        x_all.append(k)
                        y_all.append(i)
                        z_all.append(j)
                        u_all.append(u)
                        v_all.append(v)
                        w_all.append(w)
                
                else:
                    x_all.append(k)
                    y_all.append(i)
                    z_all.append(j)
                    u_all.append(u)
                    v_all.append(v)
                    w_all.append(w)
                
    x_all = np.array(x_all) / a
    y_all = np.array(y_all) / a
    z_all = np.array(z_all) / a
    u_all = np.array(u_all)
    v_all = np.array(v_all)
    w_all = np.array(w_all)
                
    return x_all, y_all, z_all, u_all, v_all, w_all

def B_rot_fun(pos, a, g, h, n, R):
    
    xyz = np.array([pos[0], pos[1], pos[2]])
    
    x_r, y_r, z_r = np.matmul(np.linalg.inv(R), xyz)
    r, theta, phi = Cart_to_Sph(x_r, y_r, z_r)
    
    frac = a/r
    
    B_r = 0
    B_theta = 0
    B_phi = 0
    
    for l in range(1, n+1):
        frac2 = frac ** (l+2)
        
        for m in range(l+1):
            
            B_r += (l+1) * frac2 * ((g[l-1][m] * np.cos(m * phi)) \
                + (h[l-1][m] * np.sin(m * phi))) * Pnm(l, m, theta)
            
            B_theta -= frac2 * ((g[l-1][m] * np.cos(m * phi)) \
                    + (h[l-1][m] * np.sin(m * phi))) * dPnm(l, m, theta)
            
            B_phi += (1 / np.sin(theta)) * m * frac2 * ((g[l-1][m] * np.sin(m * phi)) \
                    - (h[l-1][m] * np.cos(m * phi))) * Pnm(l, m, theta)
    
    B = np.array([B_r, B_theta, B_phi])        
        
    u = (np.sin(theta)*np.cos(phi) * B[0]) + (np.cos(theta)*np.cos(phi) * B[1])\
        + (-np.sin(phi) * B[2])
    v = (np.sin(theta)*np.sin(phi) * B[0]) + (np.cos(theta)*np.sin(phi) * B[1])\
        +  (np.cos(phi) * B[2])
    w = (np.cos(theta) * B[0]) + (-np.sin(theta) * B[1])
    
    B_rot = np.array([u, v, w])
    B_f = np.matmul(R, B_rot)
                
    return B_f

def B_aligned_cart(x, z, a, args, n, R, phi_rot):
   
    x_all = []
    y_all = []
    z_all = []
    u_all = []
    v_all = []
    w_all = []
    
    if n == 1:
    
        for k in x:
            for j in z:
                xyz = np.array([k, 0, j])
                R_z = np.array([[np.cos(-phi_rot), - np.sin(-phi_rot), 0],
                               [np.sin(-phi_rot), np.cos(-phi_rot), 0],
                               [0, 0, 1]])
                
                xyz2 = np.matmul(R_z, xyz)
                
                x_r, y_r, z_r = np.matmul(np.linalg.inv(R), xyz2)
                r, theta, phi = Cart_to_Sph(x_r, y_r, z_r)
                
                r01 = (2 * ((a / r) ** 3)) * (args[0] * Pnm(1, 0, theta))
                r11 = (2 * ((a / r) ** 3)) * ((args[1] * np.cos(phi)) + \
                                    (args[2] * np.sin(phi))) * Pnm(1, 1, theta)
                
                B_r = r01 + r11
                
                theta01 = ((a / r) ** 3) * (args[0] * (dPnm(1, 0, theta)))
                theta11 = ((a / r) ** 3) * ((args[1] * np.cos(phi)) + \
                                  (args[2] * np.sin(phi))) * dPnm(1, 1, theta)
                
                B_theta = - (theta01 + theta11)
                
                phi01 = 0
                phi11 = ((a / r) ** (3)) * ((args[1] * np.sin(phi)) - (args[2] * np.cos(phi))) \
                                                * Pnm(1, 1, theta)
                
                B_phi = (1 / np.sin(theta)) * (phi01 + phi11)
                
                B = np.array([B_r, B_theta, B_phi])        
        
                u = (np.sin(theta)*np.cos(phi) * B[0]) + (np.cos(theta)*np.cos(phi) * B[1])\
                    + (-np.sin(phi) * B[2])
                v = (np.sin(theta)*np.sin(phi) * B[0]) + (np.cos(theta)*np.sin(phi) * B[1])\
                    +  (np.cos(phi) * B[2])
                w = (np.cos(theta) * B[0]) + (-np.sin(theta) * B[1])   
                
                B_rot = np.array([u, v, w])
                B_f = np.matmul(R, B_rot)
                
                B_f = np.matmul(np.linalg.inv(R_z), B_f)
                
                if np.sqrt((k * k) + (j * j)) > a:
                    
                    x_all.append(k / a)
                    y_all.append(xyz2[1] / a)
                    z_all.append(j / a)
                    u_all.append(B_f[0])
                    v_all.append(B_f[1])
                    w_all.append(B_f[2])
                
                else:
                    
                    x_all.append(float("NaN"))
                    y_all.append(float("NaN"))
                    z_all.append(float("NaN"))
                    u_all.append(float("NaN"))
                    v_all.append(float("NaN"))
                    w_all.append(float("NaN"))
                    
    if n == 2:
        for k in x:
            for j in z:
                xyz = np.array([k, 0, j])
                R_z = np.array([[np.cos(phi_rot), - np.sin(phi_rot), 0],
                               [np.sin(phi_rot), np.cos(phi_rot), 0],
                               [0, 0, 1]])
                
                xyz2 = np.matmul(R_z, xyz)
                
                x_r, y_r, z_r = np.matmul(np.linalg.inv(R), xyz2)
                r, theta, phi = Cart_to_Sph(x_r, y_r, z_r)
                
                r01 = (2 * ((a / r) ** 3)) * (args[0] * Pnm(1, 0, theta))
                r11 = (2 * ((a / r) ** 3)) * ((args[1] * np.cos(phi)) + \
                                    (args[2] * np.sin(phi))) * Pnm(1, 1, theta)
                r02 = (3 * ((a / r) ** 4)) * (args[3] * Pnm(2, 0, theta))
                r12 = (3 * ((a / r) ** 4)) * ((args[4] * np.cos(phi)) + \
                                    (args[5] * np.sin(phi))) * Pnm(2, 1, theta)
                r22 = (3 * ((a / r) ** 4)) * ((args[6] * np.cos(2 * phi)) + \
                                    (args[7] * np.sin(2 * phi))) * Pnm(2, 2, theta)
                
                B_r = r01 + r11 + r02 + r12 + r22
                
                theta01 = ((a / r) ** 3) * (args[0] * (dPnm(1, 0, theta)))
                theta11 = ((a / r) ** 3) * ((args[1] * np.cos(phi)) + \
                                  (args[2] * np.sin(phi))) * dPnm(1, 1, theta)
                theta02 = ((a / r) ** 4) * (args[3] * dPnm(2, 0, theta))
                theta12 = ((a / r) ** 4) * ((args[4] * np.cos(phi)) + \
                                    (args[5] * np.sin(phi))) * dPnm(2, 1, theta)
                theta22 = ((a / r) ** 4) * ((args[6] * np.cos(2 * phi)) + \
                                    (args[7] * np.sin(2 * phi))) * dPnm(2, 2, theta)
                
                B_theta = - (theta01 + theta11 + theta02 + theta12 + theta22)
                
                phi01 = 0
                phi11 = ((a / r) ** (3)) * ((args[1] * np.sin(phi)) - (args[2] * np.cos(phi))) \
                                                * Pnm(1, 1, theta)
                phi02 = 0
                phi12 = ((a / r) ** 4) * ((args[4] * np.sin(phi)) - \
                                    (args[5] * np.cos(phi))) * Pnm(2, 1, theta)
                phi22 = (2 * ((a / r) ** 4)) * ((args[6] * np.sin(2 * phi)) - \
                                    (args[7] * np.cos(2 * phi))) * Pnm(2, 2, theta)
                
                B_phi = (1 / np.sin(theta)) * (phi01 + phi11 + phi02 + phi12 + phi22)
                
                B = np.array([B_r, B_theta, B_phi])       
                
                u = (np.sin(theta)*np.cos(phi) * B[0]) + (np.cos(theta)*np.cos(phi) * B[1])\
                    + (-np.sin(phi) * B[2])
                v = (np.sin(theta)*np.sin(phi) * B[0]) + (np.cos(theta)*np.sin(phi) * B[1])\
                    +  (np.cos(phi) * B[2])
                w = (np.cos(theta) * B[0]) + (-np.sin(theta) * B[1])
                
                B_rot = np.array([u, v, w])
                B_f = np.matmul(R, B_rot)
                
                B_f = np.matmul(np.linalg.inv(R_z), B_f)
                
                if np.sqrt((k * k) + (j * j)) > a:
                    
                    x_all.append(k / a)
                    y_all.append(xyz2[1] / a)
                    z_all.append(j / a)
                    u_all.append(B_f[0])
                    v_all.append(B_f[1])
                    w_all.append(B_f[2])
                    
                else:
                    
                    x_all.append(float("NaN"))
                    y_all.append(float("NaN"))
                    z_all.append(float("NaN"))
                    u_all.append(float("NaN"))
                    v_all.append(float("NaN"))
                    w_all.append(float("NaN"))
                    
    x_all = np.array(x_all).reshape((len(x), len(z))).T
    y_all = np.array(y_all).reshape((len(x), len(z))).T
    z_all = np.array(z_all).reshape((len(x), len(z))).T
    u_all = np.array(u_all).reshape((len(x), len(z))).T
    v_all = np.array(v_all).reshape((len(x), len(z))).T
    w_all = np.array(w_all).reshape((len(x), len(z))).T
    
    return x_all, y_all, z_all, u_all, v_all, w_all

def B_mag_cart(x1, z1, a, args, n, R, phi, matrix = True):
    
    x, y, z, u, v, w = B_aligned_cart(x1, z1, a, args, n, R, phi)
    
    B_mag = np.sqrt((u * u) + (v * v) + (w * w))

    x = x.reshape((len(x1), len(z1))).T
    y = y.reshape((len(x1), len(z1))).T
    z = z.reshape((len(x1), len(z1))).T
    u = u.reshape((len(x1), len(z1))).T
    v = v.reshape((len(x1), len(z1))).T
    w = w.reshape((len(x1), len(z1))).T
    if matrix == True:
        B_mag = B_mag.reshape((len(x1), len(z1))).T
    else:
        B_mag = B_mag.T
    return B_mag

def B_spin_aligned(r, theta, phi, a, args, n, R): #should be field aligned?
    x, y, z, u, v, w = Get_B_sph_rot(r, theta, phi, a, args, n, R)

    R_z = np.array([[np.cos(-phi[0]), - np.sin(-phi[0]), 0],
                   [np.sin(-phi[0]), np.cos(-phi[0]), 0],
                   [0, 0, 1]])

    xyz = np.array([x, y, z])
    x, y, z = np.matmul(R_z, xyz)
    uvw = np.array([u, v, w])
    u, v, w = np.matmul(R_z, uvw)
    
    return x, y, z, u, v, w

def B_magnitude(r, theta, phi, a, args, n, R):
    x, y, z, u, v, w = B_spin_aligned(r, theta, phi, a, args, n, R)

    B_mag = np.sqrt((u * u) + (v * v) + (w * w))
    
    xyz = np.array([x, y, z])
    uvw = np.array([u, v, w])
    
    return B_mag, xyz, uvw 

def Get_maxB_ratio(r, theta, phi, a, args, R, plot):
    B_all = []
    
    for i in phi:
        
        uq, vq, wq = B_magnitude(r, theta, [i], a, args, 2, R)[2]
        ud, vd, wd = B_magnitude(r, theta, [i], a, args, 1, R)[2] 
    
        uq_d = uq - ud
        vq_d = vq - vd
        wq_d = wq - wd
        top = np.sqrt((uq_d * uq_d) + (vq_d * vq_d) + (wq_d * wq_d))
        bottom = np.sqrt((ud * ud) + (vd * vd) + (wd * wd))

        B_ratio = abs(top) / abs(bottom)
        
        B_all.append(max(B_ratio))
    
    if plot == True:
        plt.plot(phi * 180 / np.pi, B_all)
        plt.xlabel('Longitude (degrees)')
        plt.ylabel('Max ratio')
        plt.show()    
    
    return B_all

def B_Lshell(L, theta_in, phi_in, a, args, n, R):
  
    x_all = []
    y_all = []
    z_all = []
    u_all = []
    v_all = []
    w_all = []
    
    
    if n == 1:
        for j in theta_in:
            for k in phi_in:
                
                r_L = L * np.sin(j) * np.sin(j) * a
                
                x, y, z = Sph_to_Cart(r_L, j, k)
                xyz = np.array([x, y, z])
                x_r, y_r, z_r = np.matmul(np.linalg.inv(R), xyz)
                r, theta, phi = Cart_to_Sph(x_r, y_r, z_r)
                
                r01 = (2 * ((a / r) ** 3)) * (args[0] * Pnm(1, 0, theta))
                r11 = (2 * ((a / r) ** 3)) * ((args[1] * np.cos(phi)) + \
                                    (args[2] * np.sin(phi))) * Pnm(1, 1, theta)
                
                B_r = r01 + r11
                
                theta01 = ((a / r) ** 3) * (args[0] * (dPnm(1, 0, theta)))
                theta11 = ((a / r) ** 3) * ((args[1] * np.cos(phi)) + \
                                  (args[2] * np.sin(phi))) * dPnm(1, 1, theta)
                
                B_theta = - (theta01 + theta11)
                
                phi01 = 0
                phi11 = ((a / r) ** (3)) * ((args[1] * np.sin(phi)) - (args[2] * np.cos(phi))) \
                                                * Pnm(1, 1, theta)
                
                B_phi = (1 / np.sin(theta)) * (phi01 + phi11)
                
                B = np.array([B_r, B_theta, B_phi])        
        
                u = (np.sin(theta)*np.cos(phi) * B[0]) + (np.cos(theta)*np.cos(phi) * B[1])\
                    + (-np.sin(phi) * B[2])
                v = (np.sin(theta)*np.sin(phi) * B[0]) + (np.cos(theta)*np.sin(phi) * B[1])\
                    +  (np.cos(phi) * B[2])
                w = (np.cos(theta) * B[0]) + (-np.sin(theta) * B[1])   
                
                B_rot = np.array([u, v, w])
                B_f = np.matmul(R, B_rot)
                
                x_all.append(x)
                y_all.append(y)
                z_all.append(z)
                u_all.append(B_f[0])
                v_all.append(B_f[1])
                w_all.append(B_f[2])
    
    if n == 2:
        for j in theta_in:
            for k in phi_in:
                
                r_L = L * np.sin(j) * np.sin(j) * a
                
                x, y, z = Sph_to_Cart(r_L, j, k)
                xyz = np.array([x, y, z])
                x_r, y_r, z_r = np.matmul(np.linalg.inv(R), xyz)
                r, theta, phi = Cart_to_Sph(x_r, y_r, z_r)
                
                r01 = (2 * ((a / r) ** 3)) * (args[0] * Pnm(1, 0, theta))
                r11 = (2 * ((a / r) ** 3)) * ((args[1] * np.cos(phi)) + \
                                    (args[2] * np.sin(phi))) * Pnm(1, 1, theta)
                r02 = (3 * ((a / r) ** 4)) * (args[3] * Pnm(2, 0, theta))
                r12 = (3 * ((a / r) ** 4)) * ((args[4] * np.cos(phi)) + \
                                    (args[5] * np.sin(phi))) * Pnm(2, 1, theta)
                r22 = (3 * ((a / r) ** 4)) * ((args[6] * np.cos(2 * phi)) + \
                                    (args[7] * np.sin(2 * phi))) * Pnm(2, 2, theta)
                
                B_r = r01 + r11 + r02 + r12 + r22
                
                theta01 = ((a / r) ** 3) * (args[0] * (dPnm(1, 0, theta)))
                theta11 = ((a / r) ** 3) * ((args[1] * np.cos(phi)) + \
                                  (args[2] * np.sin(phi))) * dPnm(1, 1, theta)
                theta02 = ((a / r) ** 4) * (args[3] * dPnm(2, 0, theta))
                theta12 = ((a / r) ** 4) * ((args[4] * np.cos(phi)) + \
                                    (args[5] * np.sin(phi))) * dPnm(2, 1, theta)
                theta22 = ((a / r) ** 4) * ((args[6] * np.cos(2 * phi)) + \
                                    (args[7] * np.sin(2 * phi))) * dPnm(2, 2, theta)
                
                B_theta = - (theta01 + theta11 + theta02 + theta12 + theta22)
                
                phi01 = 0
                phi11 = ((a / r) ** (3)) * ((args[1] * np.sin(phi)) - (args[2] * np.cos(phi))) \
                                                * Pnm(1, 1, theta)
                phi02 = 0
                phi12 = ((a / r) ** 4) * ((args[4] * np.sin(phi)) - \
                                    (args[5] * np.cos(phi))) * Pnm(2, 1, theta)
                phi22 = (2 * ((a / r) ** 4)) * ((args[6] * np.sin(2 * phi)) - \
                                    (args[7] * np.cos(2 * phi))) * Pnm(2, 2, theta)
                
                B_phi = (1 / np.sin(theta)) * (phi01 + phi11 + phi02 + phi12 + phi22)
                
                B = np.array([B_r, B_theta, B_phi])       
                
                u = (np.sin(theta)*np.cos(phi) * B[0]) + (np.cos(theta)*np.cos(phi) * B[1])\
                    + (-np.sin(phi) * B[2])
                v = (np.sin(theta)*np.sin(phi) * B[0]) + (np.cos(theta)*np.sin(phi) * B[1])\
                    +  (np.cos(phi) * B[2])
                w = (np.cos(theta) * B[0]) + (-np.sin(theta) * B[1])
                
                B_rot = np.array([u, v, w])
                B_f = np.matmul(R, B_rot)
                
                x_all.append(x)
                y_all.append(y)
                z_all.append(z)
                u_all.append(B_f[0])
                v_all.append(B_f[1])
                w_all.append(B_f[2])
                    
    x_all = np.array(x_all) / a
    y_all = np.array(y_all) / a
    z_all = np.array(z_all) / a
    u_all = np.array(u_all)
    v_all = np.array(v_all)
    w_all = np.array(w_all)
                    
    return x_all, y_all, z_all, u_all, v_all, w_all

def gradB_fun(x, y, z, a, g, h, n, R):
    
    x_all = []
    y_all = []
    z_all = []
    du_all = []
    dv_all = []
    dw_all = []
    u_all = []
    v_all = []
    w_all = []
    
    for k in x:
        for i in y:
            for j in z:
                xyz = np.array([k, i, j])
                x_r, y_r, z_r = np.matmul(np.linalg.inv(R), xyz)
                r, theta, phi = Cart_to_Sph(x_r, y_r, z_r)
                
                frac = a/r
                
                Brdr = 0
                Brdtheta = 0
                Brdphi = 0
                Bthetadr = 0
                Bthetadtheta = 0
                Bthetadphi = 0
                Bphidr = 0
                Bphidtheta = 0
                Bphidphi = 0
                
                for l in range(1, n+1):
                    frac2 = frac ** (l+2)
                    
                    for m in range(l+1):
                        
                        cos = np.cos(m * phi)
                        sin = np.sin(m * phi)
                        
                        Brdr -= (l+1) * (l+2) * (a ** (l+2)) * (r ** (-l-3)) * ((g[l-1][m] * cos) \
                            + (h[l-1][m] * sin)) * Pnm(l, m, theta)
                        
                        Brdtheta += (l+1) * frac2 * ((g[l-1][m] * cos) \
                            + (h[l-1][m] * sin)) * dPnm(l, m, theta)
                            
                        Brdphi += (l+1) * frac2 * ((g[l-1][m] * (-m) * sin) \
                            + (h[l-1][m] * m * cos)) * Pnm(l, m, theta)
                        
                        Bthetadr += (l+2) * (a ** (l+2)) * (r ** (-l-3)) * ((g[l-1][m] * cos) \
                                + (h[l-1][m] * sin)) * dPnm(l, m, theta)
                        
                        Bthetadtheta -= frac2 * ((g[l-1][m] * cos) \
                                + (h[l-1][m] * sin)) * d2Pnm(l, m, theta)
                            
                        Bthetadphi -= frac2 * ((g[l-1][m] * (-m) * sin) \
                                + (h[l-1][m] * m * cos)) * dPnm(l, m, theta)
                            
                        Bphidr -= (1 / np.sin(theta)) * m * (l+2) * (a ** (l+2)) * (r ** (-l-3))\
                            * ((g[l-1][m] * sin) - (h[l-1][m] * cos)) * Pnm(l, m, theta)
                        
                        Bphidtheta += m * frac2 * ((g[l-1][m] * sin) - (h[l-1][m] * cos)) * \
                            (((np.sin(theta) * dPnm(l, m, theta)) - (Pnm(l, m, theta) * \
                            np.cos(theta))) / (np.sin(theta) ** 2))
                        
                        Bphidphi += (1 / np.sin(theta)) * m * frac2 * ((g[l-1][m] * m * cos) \
                                + (h[l-1][m] * m * sin)) * Pnm(l, m, theta)
                
                B = B_sph(r, theta, phi, a, g, h, n)
                
                dBdr = ((B[0] * Brdr) + (B[1] * Bthetadr) + (B[2] * Bphidr)) / np.linalg.norm(B)           
                dBdtheta = ((B[0] * Brdtheta) + (B[1] * Bthetadtheta) + (B[2] * Bphidtheta)) / np.linalg.norm(B) 
                dBdphi = ((B[0] * Brdphi) + (B[1] * Bthetadphi) + (B[2] * Bphidphi)) / np.linalg.norm(B)
                
                dB = np.array([dBdr, dBdtheta / r, dBdphi / (r * np.sin(theta))])        
    
                du = (np.sin(theta)*np.cos(phi) * dB[0]) + (np.cos(theta)*np.cos(phi) * dB[1])\
                    + (-np.sin(phi) * dB[2])
                dv = (np.sin(theta)*np.sin(phi) * dB[0]) + (np.cos(theta)*np.sin(phi) * dB[1])\
                    +  (np.cos(phi) * dB[2])
                dw = (np.cos(theta) * dB[0]) + (-np.sin(theta) * dB[1])
                
                dB_rot = np.array([du, dv, dw])
                dB_f = np.matmul(R, dB_rot)
                
                u = (np.sin(theta)*np.cos(phi) * B[0]) + (np.cos(theta)*np.cos(phi) * B[1])\
                    + (-np.sin(phi) * B[2])
                v = (np.sin(theta)*np.sin(phi) * B[0]) + (np.cos(theta)*np.sin(phi) * B[1])\
                    +  (np.cos(phi) * B[2])
                w = (np.cos(theta) * B[0]) + (-np.sin(theta) * B[1])
                
                B_rot = np.array([u, v, w])
                B_f = np.matmul(R, B_rot)
                
                if np.sqrt((k * k) + (i * i) + (j * j)) > a:
                    
                    x_all.append(k / a)
                    y_all.append(i / a)
                    z_all.append(j / a)
                    du_all.append(dB_f[0])
                    dv_all.append(dB_f[1])
                    dw_all.append(dB_f[2])
                    u_all.append(B_f[0])
                    v_all.append(B_f[1])
                    w_all.append(B_f[2])
                
                else:
                    x_all.append(float("NaN"))
                    y_all.append(float("NaN"))
                    z_all.append(float("NaN"))
                    du_all.append(float("NaN"))
                    dv_all.append(float("NaN"))
                    dw_all.append(float("NaN"))
                    u_all.append(float("NaN"))
                    v_all.append(float("NaN"))
                    w_all.append(float("NaN"))
    
    x_all = np.array(x_all)
    y_all = np.array(y_all)
    z_all = np.array(z_all)
    du_all = np.array(du_all)
    dv_all = np.array(dv_all)
    dw_all = np.array(dw_all)
    u_all = np.array(u_all)
    v_all = np.array(v_all)
    w_all = np.array(w_all)
    
    coords = np.array([x_all, y_all, z_all])
    dB = np.array([du_all, dv_all, dw_all])
    B = np.array([u_all, v_all, w_all])
    
    return coords, dB, B

def gradB_long(x, z, a, g, h, n, R, phi_rot):
    
    x_all = []
    y_all = []
    z_all = []
    du_all = []
    dv_all = []
    dw_all = []
    u_all = []
    v_all = []
    w_all = []
    
    for k in x:
        for j in z:
            xyz = np.array([k, 0, j])
            R_z = np.array([[np.cos(phi_rot), - np.sin(phi_rot), 0],
                           [np.sin(phi_rot), np.cos(phi_rot), 0],
                           [0, 0, 1]])
            
            xyz2 = np.matmul(R_z, xyz)
            
            x_r, y_r, z_r = np.matmul(np.linalg.inv(R), xyz2)
            r, theta, phi = Cart_to_Sph(x_r, y_r, z_r)
            
            frac = a/r
            
            Brdr = 0
            Brdtheta = 0
            Brdphi = 0
            Bthetadr = 0
            Bthetadtheta = 0
            Bthetadphi = 0
            Bphidr = 0
            Bphidtheta = 0
            Bphidphi = 0
            
            for l in range(1, n+1):
                frac2 = frac ** (l+2)
                
                for m in range(l+1):
                    
                    cos = np.cos(m * phi)
                    sin = np.sin(m * phi)
                    
                    Brdr -= (l+1) * (l+2) * (a ** (l+2)) * (r ** (-l-3)) * ((g[l-1][m] * cos) \
                        + (h[l-1][m] * sin)) * Pnm(l, m, theta)
                    
                    Brdtheta += (l+1) * frac2 * ((g[l-1][m] * cos) \
                        + (h[l-1][m] * sin)) * dPnm(l, m, theta)
                        
                    Brdphi += (l+1) * frac2 * ((g[l-1][m] * (-m) * sin) \
                        + (h[l-1][m] * m * cos)) * Pnm(l, m, theta)
                    
                    Bthetadr += (l+2) * (a ** (l+2)) * (r ** (-l-3)) * ((g[l-1][m] * cos) \
                            + (h[l-1][m] * sin)) * dPnm(l, m, theta)
                    
                    Bthetadtheta -= frac2 * ((g[l-1][m] * cos) \
                            + (h[l-1][m] * sin)) * d2Pnm(l, m, theta)
                        
                    Bthetadphi -= frac2 * ((g[l-1][m] * (-m) * sin) \
                            + (h[l-1][m] * m * cos)) * dPnm(l, m, theta)
                        
                    Bphidr -= (1 / np.sin(theta)) * m * (l+2) * (a ** (l+2)) * (r ** (-l-3))\
                        * ((g[l-1][m] * sin) - (h[l-1][m] * cos)) * Pnm(l, m, theta)
                    
                    Bphidtheta += m * frac2 * ((g[l-1][m] * sin) - (h[l-1][m] * cos)) * \
                        (((np.sin(theta) * dPnm(l, m, theta)) - (Pnm(l, m, theta) * \
                        np.cos(theta))) / (np.sin(theta) ** 2))
                    
                    Bphidphi += (1 / np.sin(theta)) * m * frac2 * ((g[l-1][m] * m * cos) \
                            + (h[l-1][m] * m * sin)) * Pnm(l, m, theta)
            
            B = B_sph(r, theta, phi, a, g, h, n)
            
            dBdr = ((B[0] * Brdr) + (B[1] * Bthetadr) + (B[2] * Bphidr)) / np.linalg.norm(B)           
            dBdtheta = ((B[0] * Brdtheta) + (B[1] * Bthetadtheta) + (B[2] * Bphidtheta)) / np.linalg.norm(B) 
            dBdphi = ((B[0] * Brdphi) + (B[1] * Bthetadphi) + (B[2] * Bphidphi)) / np.linalg.norm(B)
            
            dB = np.array([dBdr, dBdtheta / r, dBdphi / (r * np.sin(theta))])        

            du = (np.sin(theta)*np.cos(phi) * dB[0]) + (np.cos(theta)*np.cos(phi) * dB[1])\
                + (-np.sin(phi) * dB[2])
            dv = (np.sin(theta)*np.sin(phi) * dB[0]) + (np.cos(theta)*np.sin(phi) * dB[1])\
                +  (np.cos(phi) * dB[2])
            dw = (np.cos(theta) * dB[0]) + (-np.sin(theta) * dB[1])
            
            dB_rot = np.array([du, dv, dw])
            dB_f = np.matmul(R, dB_rot)
            
            dB_f = np.matmul(np.linalg.inv(R_z), dB_f)
            
            u = (np.sin(theta)*np.cos(phi) * B[0]) + (np.cos(theta)*np.cos(phi) * B[1])\
                + (-np.sin(phi) * B[2])
            v = (np.sin(theta)*np.sin(phi) * B[0]) + (np.cos(theta)*np.sin(phi) * B[1])\
                +  (np.cos(phi) * B[2])
            w = (np.cos(theta) * B[0]) + (-np.sin(theta) * B[1])
            
            B_rot = np.array([u, v, w])
            B_f = np.matmul(R, B_rot)
            
            B_f = np.matmul(np.linalg.inv(R_z), B_f)
            
            if np.sqrt((k * k) + (j * j)) > a:
                
                x_all.append(k / a)
                y_all.append(xyz2[1] / a)
                z_all.append(j / a)
                du_all.append(dB_f[0])
                dv_all.append(dB_f[1])
                dw_all.append(dB_f[2])
                u_all.append(B_f[0])
                v_all.append(B_f[1])
                w_all.append(B_f[2])
            
            else:
                x_all.append(float("NaN"))
                y_all.append(float("NaN"))
                z_all.append(float("NaN"))
                du_all.append(float("NaN"))
                dv_all.append(float("NaN"))
                dw_all.append(float("NaN"))
                u_all.append(float("NaN"))
                v_all.append(float("NaN"))
                w_all.append(float("NaN"))
    
    x_all = np.array(x_all).reshape((len(x), len(z))).T
    y_all = np.array(y_all).reshape((len(x), len(z))).T
    z_all = np.array(z_all).reshape((len(x), len(z))).T
    du_all = np.array(du_all).reshape((len(x), len(z))).T
    dv_all = np.array(dv_all).reshape((len(x), len(z))).T
    dw_all = np.array(dw_all).reshape((len(x), len(z))).T
    u_all = np.array(u_all).reshape((len(x), len(z))).T
    v_all = np.array(v_all).reshape((len(x), len(z))).T
    w_all = np.array(w_all).reshape((len(x), len(z))).T
    
    coords = np.array([x_all, y_all, z_all])
    dB = np.array([du_all, dv_all, dw_all])
    B = np.array([u_all, v_all, w_all])
    
    return coords, dB, B

def dB_Lshell(L, theta_in, phi_in, a, g, h, n, R):
  
    x_all = []
    y_all = []
    z_all = []
    du_all = []
    dv_all = []
    dw_all = []
    u_all = []
    v_all = []
    w_all = []
    normal_all = []
    phi_hat_all = []
    sigma_all = []
    
    
    for j in theta_in:
        for k in phi_in:
            
            r_L = L * np.sin(j) * np.sin(j) * a
            
            x, y, z = Sph_to_Cart(r_L, j, k)
            xyz = np.array([x, y, z])
            x_r, y_r, z_r = np.matmul(np.linalg.inv(R), xyz)
            r, theta, phi = Cart_to_Sph(x_r, y_r, z_r)
            
            frac = a/r
            
            Brdr = 0
            Brdtheta = 0
            Brdphi = 0
            Bthetadr = 0
            Bthetadtheta = 0
            Bthetadphi = 0
            Bphidr = 0
            Bphidtheta = 0
            Bphidphi = 0
            
            for l in range(1, n+1):
                frac2 = frac ** (l+2)
                
                for m in range(l+1):
                    
                    cos = np.cos(m * phi)
                    sin = np.sin(m * phi)
                    
                    Brdr -= (l+1) * (l+2) * (a ** (l+2)) * (r ** (-l-3)) * ((g[l-1][m] * cos) \
                        + (h[l-1][m] * sin)) * Pnm(l, m, theta)
                    
                    Brdtheta += (l+1) * frac2 * ((g[l-1][m] * cos) \
                        + (h[l-1][m] * sin)) * dPnm(l, m, theta)
                        
                    Brdphi += (l+1) * frac2 * ((g[l-1][m] * (-m) * sin) \
                        + (h[l-1][m] * m * cos)) * Pnm(l, m, theta)
                    
                    Bthetadr += (l+2) * (a ** (l+2)) * (r ** (-l-3)) * ((g[l-1][m] * cos) \
                            + (h[l-1][m] * sin)) * dPnm(l, m, theta)
                    
                    Bthetadtheta -= frac2 * ((g[l-1][m] * cos) \
                            + (h[l-1][m] * sin)) * d2Pnm(l, m, theta)
                        
                    Bthetadphi -= frac2 * ((g[l-1][m] * (-m) * sin) \
                            + (h[l-1][m] * m * cos)) * dPnm(l, m, theta)
                        
                    Bphidr -= (1 / np.sin(theta)) * m * (l+2) * (a ** (l+2)) * (r ** (-l-3))\
                        * ((g[l-1][m] * sin) - (h[l-1][m] * cos)) * Pnm(l, m, theta)
                    
                    Bphidtheta += m * frac2 * ((g[l-1][m] * sin) - (h[l-1][m] * cos)) * \
                        (((np.sin(theta) * dPnm(l, m, theta)) - (Pnm(l, m, theta) * \
                        np.cos(theta))) / (np.sin(theta) ** 2))
                    
                    Bphidphi += (1 / np.sin(theta)) * m * frac2 * ((g[l-1][m] * m * cos) \
                            + (h[l-1][m] * m * sin)) * Pnm(l, m, theta)
            
            B = B_sph(r, theta, phi, a, g, h, n)
            
            dBdr = ((B[0] * Brdr) + (B[1] * Bthetadr) + (B[2] * Bphidr)) / np.linalg.norm(B)           
            dBdtheta = ((B[0] * Brdtheta) + (B[1] * Bthetadtheta) + (B[2] * Bphidtheta)) / np.linalg.norm(B) 
            dBdphi = ((B[0] * Brdphi) + (B[1] * Bthetadphi) + (B[2] * Bphidphi)) / np.linalg.norm(B)
            
            dB = np.array([dBdr, dBdtheta / r, dBdphi / (r * np.sin(theta))])        

            du = (np.sin(theta)*np.cos(phi) * dB[0]) + (np.cos(theta)*np.cos(phi) * dB[1])\
                + (-np.sin(phi) * dB[2])
            dv = (np.sin(theta)*np.sin(phi) * dB[0]) + (np.cos(theta)*np.sin(phi) * dB[1])\
                +  (np.cos(phi) * dB[2])
            dw = (np.cos(theta) * dB[0]) + (-np.sin(theta) * dB[1])
            
            dB_rot = np.array([du, dv, dw])
            dB_f = np.matmul(R, dB_rot)
            
            u = (np.sin(theta)*np.cos(phi) * B[0]) + (np.cos(theta)*np.cos(phi) * B[1])\
                + (-np.sin(phi) * B[2])
            v = (np.sin(theta)*np.sin(phi) * B[0]) + (np.cos(theta)*np.sin(phi) * B[1])\
                +  (np.cos(phi) * B[2])
            w = (np.cos(theta) * B[0]) + (-np.sin(theta) * B[1])
            
            B_rot = np.array([u, v, w])
            B_f = np.matmul(R, B_rot)
            
            xn =  x - ((2 * z * z * x) / ((x * x) + (y * y)))
            yn =  y - ((2 * z * z * y) / ((x * x) + (y * y)))
            zn = z + (2 * z)

            normal = np.array([xn, yn, zn])
            normal = normal / np.linalg.norm(normal)
            
            phi_hat = np.array([-y, x, 0]) / np.sqrt((x * x) + (y * y))
            
            sigma = np.cross(phi_hat, normal)
            sigma = sigma / np.linalg.norm(sigma)
            
            x_all.append(x)
            y_all.append(y)
            z_all.append(z)
            du_all.append(dB_f[0])
            dv_all.append(dB_f[1])
            dw_all.append(dB_f[2])
            u_all.append(B_f[0])
            v_all.append(B_f[1])
            w_all.append(B_f[2])
            normal_all.append(normal)
            phi_hat_all.append(phi_hat)
            sigma_all.append(sigma)
                    
    x_all = np.array(x_all) / a
    y_all = np.array(y_all) / a
    z_all = np.array(z_all) / a
    du_all = np.array(du_all)
    dv_all = np.array(dv_all)
    dw_all = np.array(dw_all)
    u_all = np.array(u_all)
    v_all = np.array(v_all)
    w_all = np.array(w_all)
                    
    coords = np.array([x_all, y_all, z_all])
    dB = np.array([du_all, dv_all, dw_all])
    B = np.array([u_all, v_all, w_all])
    
    return coords, dB, B, normal_all, phi_hat_all, sigma_all

def TwoD_plot(x1, x2, y1, y2, plane):
    """
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
    """
    if plane == 'x':
        y = x1
        z = x2
        v = y1
        w = y2
        
        
        Circle1 = plt.Circle((0, 0), 1, color = 'red', fill = False)
        fig, ax = plt.subplots()
        ax.add_patch(Circle1)
        ax.quiver(y, z, v, w, pivot = 'mid')
        plt.title('at $x/a$ = 0')
        plt.xlabel('y/a')
        plt.ylabel('z/a')
        #plt.xlim(-2, 2)
        #plt.ylim(-2, 2)
        #plt.savefig('Rot_dipole/z=0')
        plt.show()
        
    elif plane == 'y':
        x = x1
        z = x2
        u = y1
        w = y2
        
        Circle1 = plt.Circle((0, 0), 1, color = 'red', fill = False)
        fig, ax = plt.subplots()
        ax.add_patch(Circle1)
        ax.quiver(x, z, u, w, pivot = 'mid')
        plt.title('at $y/a$ = 0')
        plt.xlabel('x/a')
        plt.ylabel('z/a')
        #plt.xlim(-2, 2)
        #plt.ylim(-2, 2)
        #plt.savefig('Rot_dipole/z=0')
        plt.show()
        
    elif plane == 'z':
        x = x1
        y = x2
        u = y1
        v = y2
        
        Circle1 = plt.Circle((0, 0), 1, color = 'red', fill = False)
        fig, ax = plt.subplots()
        ax.add_patch(Circle1)
        ax.quiver(x, y, u, v, pivot = 'mid')
        plt.title('at $z/a$ = 0')
        plt.xlabel('x/a')
        plt.ylabel('y/a')
        #plt.xlim(-2, 2)
        #plt.ylim(-2, 2)
        #plt.savefig('Rot_dipole/z=0')
        plt.show()
    
def ThreeD_plot(x, y, z, u, v, w):
    params = {
   'axes.labelsize': 20,
   'font.size': 20,
   #'font.family': 'sans-serif', # Optionally change the font family to sans-serif
   #'font.serif': 'Arial', # Optionally change the font to Arial
   'legend.fontsize': 10,
   'xtick.labelsize': 20,
   'ytick.labelsize': 20, 
   'figure.figsize': [20, 20]
    } 
    plt.rcParams.update(params)
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    # ax.set_xlim3d(-2, 2)
    # ax.set_ylim3d(-2, 2)
    # ax.set_zlim3d(-2, 2)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    
    # zline = np.linspace(-7, 7, 10)
    # ax.plot3D(zline * 0, zline * 0, zline, 'red')
    ax.quiver(x, y, z, u, v, w, length=0.4, normalize=True)
    plt.show()

def animate_vector(num, qr, r, theta, a, args, R):
   phi = np.array([num, num + np.pi])

   x, y, z, u, v, w = Get_B_sph_rot(r, theta, phi, a, args, 2, R)
   
   R_z = np.array([[np.cos(-phi[0]), - np.sin(-phi[0]), 0],
                  [np.sin(-phi[0]), np.cos(-phi[0]), 0],
                  [0, 0, 1]])

   xyz = np.array([x, y, z])
   x, y, z = np.matmul(R_z, xyz)
   uvw = np.array([u, v, w])
   u, v, w = np.matmul(R_z, uvw)
   qr.set_UVC(u, w)
   plt.title('Phi = {} $\pi$'.format(num / np.pi))
   plt.ylabel('z/α')
   return qr,

def animate_colourmap_ratio(num, scat, r, theta, a, args, R):
   phi = [num]

   B_mag_quad, xyz_quad, uvw_quad = B_magnitude(r, theta, phi, a, args, 2, R)
   B_mag_dip, xyz_dip, uvw_dip = B_magnitude(r, theta, phi, a, args, 1, R)

   B_ratio = abs((B_mag_quad - B_mag_dip) / B_mag_dip)
   
   scat.set_array(B_ratio)
   
   plt.title('Phi = {:.1f} $\pi$'.format(num / np.pi))
   plt.ylabel('z/α')
   #plt.colorbar()
   return scat,

