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
from matplotlib import animation
import matplotlib
from tqdm import tqdm
from matplotlib import colors
#%%
#CONSTANTS    
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

args = np.array([g01, g11, h11, g02, g12, h12, g22, h22])
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

r = np.linspace(2 * a, 3 * a, 5)
theta = np.linspace(0, np.pi , 10)
phi = np.linspace(0, 2 * np.pi, 10)

x, y, z, u, v, w = Get_B_sph(r, theta, phi, a, args, 2)

ThreeD_plot(x, y, z, u, v, w)

#%%
z1 = [0]#np.linspace(-1.5 * a, 1.5 * a, 4)
y1 = np.linspace(-1.5 * a, 1.5 * a, 20)
x1 = np.linspace(-1.5 * a, 1.5 * a, 20)

x, y, z, u, v, w = Get_B_cart(x1, y1, z1, a, args, 2, True)

TwoD_plot(x, y, u, v, 'z')
            
#%%
x1 = [0]#np.linspace(-1.5 * a, 1.5 * a, 4)
y1 = np.linspace(-1.5 * a, 1.5 * a, 4)
z1 = np.linspace(-1.5 * a, 1.5 * a, 4)

g = np.array([[g01, g11], 
              [g02, g12, g22]], dtype = object)
h = np.array([[0, h11], 
              [0, h12, h22]], dtype = object)

x, y, z, u, v, w = getB_fun(x1, y1, z1, a, g, h, 2)

TwoD_plot(y, z, v, w, 'x')

#print(max(u))

#%%

a = 25600000 #Uranus radius

z1 = [0]#np.linspace(-1.5 * a, 1.5 * a, 20)
y1 = np.linspace(-3 * a, 3 * a, 20)
x1 = np.linspace(-3 * a, 3 * a, 20)

x, y, z, u, v, w = Get_B_cart(x1, y1, z1, a, args, 2, True)

TwoD_plot(x, y, u, v, 'z')

#%%
x1 = np.linspace(-2* a, 2 * a, 5)
z1 = np.linspace(-2 * a, 2 * a, 5)
y1 = np.linspace(-2 * a, 2 * a, 5)

x, y, z, u, v, w = Get_B_cart(x1, y1, z1, a, args, 2)

ThreeD_plot(x, y, z, u, v, w)

#%%

plt.rcParams["figure.figsize"] = [8, 8]
plt.rcParams["figure.autolayout"] = True

r = np.linspace(a, 2 * a, 4)
theta = np.linspace(0, np.pi, 10)
phi = np.array([0, np.pi])

x, y, z, u, v, w = B_spin_aligned(r, theta, phi, a, args, 2, R)

fig, ax = plt.subplots(1, 1)
qr = ax.quiver(x, z, u, w, color='black', pivot = 'mid')
plt.title('Phi = 0')

num_range = np.linspace(0, np.pi, 5)

anim_vector = animation.FuncAnimation(fig, animate_vector, frames = num_range, fargs=(qr, r, theta, a, args, R),
                              interval=1000, blit=False)

Circle1 = plt.Circle((0, 0), 1, color = 'red', fill = False)
ax.add_patch(Circle1)

plt.show()
#TwoD_plot(x, z, u, w, 'y')

#anim_vector.save('Harry/quadrupole_test.gif')

#%%
plt.rcParams["figure.figsize"] = [10, 10]
plt.rcParams["figure.autolayout"] = True
plt.scatter(x=x,y=z,c=B_ratio, s=40, norm=matplotlib.colors.LogNorm())
plt.colorbar()
#plt.contourf(X_mat, Z_mat, B_mag_mat)
plt.ylabel('$z / \u03B1$')
plt.quiver(x[::13], z[::13], u[::13], w[::13], scale = 1000000, pivot = 'mid')
plt.show()

#%%

r = np.linspace(2*a, 5 * a, 10)
theta = np.linspace(0, np.pi, 10)
phi = np.linspace(0, 2*np.pi, 100)#np.array([0, np.pi])

B_rall = Get_maxB_ratio(r, theta, phi, a, args, R, True) 
    
#     plt.scatter(x=x,y=z,c=B_ratio, s=150)#, norm=matplotlib.colors.LogNorm())
#     plt.colorbar()
#     plt.quiver(x[::13], z[::13], u[::13], w[::13], scale = 1000000, pivot = 'mid')
#     plt.title('Longitude = {:.1f}'.format(i  * 180 / np.pi))
    
#     plt.show()

#%%

r = np.linspace(1*a, 2 * a, 10)
theta = np.linspace(0, np.pi, 10)
phi = [0]#np.linspace(0, 2*np.pi, 50)#np.array([0, np.pi])

B_mag_quad, xyz_quad, uvw_quad = B_magnitude(r, theta, phi, a, args, 2, R)
B_mag_dip, xyz_dip, uvw_dip = B_magnitude(r, theta, phi, a, args, 1, R)

B_ratio = abs((B_mag_quad - B_mag_dip) / B_mag_dip)

fig, ax = plt.subplots(1, 1)
scat = plt.scatter(x=xyz_quad[0],y=xyz_quad[2],c=B_ratio, s=150)#, norm=matplotlib.colors.LogNorm())
#plt.title('Phi = 0')

num_range = np.linspace(0, 2 * np.pi, 11)

anim2 = animation.FuncAnimation(fig, animate_colourmap_ratio, frames = num_range, \
                                fargs=(scat, r, theta, a, args, R),interval=500, blit=False)

# Circle1 = plt.Circle((0, 0), 1, color = 'red', fill = False)
# ax.add_patch(Circle1)

plt.show()

#anim2.save('Harry/colormap_test.gif')

#%%
z1 = np.linspace(-5 * a, 5 * a, 10)
x1 = np.linspace(0, 5 * a, 10)

phi = 0#np.pi/0.5

x_quad, y_quad, z_quad, u_quad, v_quad, w_quad = B_aligned_cart(x1, z1, a, args, 2, R, phi)
x_dip, y_dip, z_dip, u_dip, v_dip, w_dip = B_aligned_cart(x1, z1, a, args, 1, R, phi)

uq_d = u_quad - u_dip
vq_d = v_quad - v_dip
wq_d = w_quad - w_dip
top = np.sqrt((uq_d * uq_d) + (vq_d * vq_d) + (wq_d * wq_d))
bottom = np.sqrt((u_dip * u_dip) + (v_dip * v_dip) + (w_dip * w_dip))

B_ratio = abs(top) / abs(bottom)

fig, ax = plt.subplots(1, 1)
mag_vec = np.sqrt((u_quad * u_quad) + (w_quad * w_quad))
u = u_quad / mag_vec
w = w_quad / mag_vec

u = np.delete(u, slice(0, 100, 2), 0)
u = np.delete(u, slice(0, 100, 2), 1)
u = np.delete(u, slice(0, 50, 2), 0)
u = np.delete(u, slice(0, 50, 2), 1)

w = np.delete(w, slice(0, 100, 2), 0)
w = np.delete(w, slice(0, 100, 2), 1)
w = np.delete(w, slice(0, 50, 2), 0)
w = np.delete(w, slice(0, 50, 2), 1)

x_quad = np.delete(x_quad, slice(0, 100, 2), 0)
x_quad = np.delete(x_quad, slice(0, 100, 2), 1)
x_quad = np.delete(x_quad, slice(0, 50, 2), 0)
x_quad = np.delete(x_quad, slice(0, 50, 2), 1)

z_quad = np.delete(z_quad, slice(0, 100, 2), 0)
z_quad = np.delete(z_quad, slice(0, 100, 2), 1)
z_quad = np.delete(z_quad, slice(0, 50, 2), 0)
z_quad = np.delete(z_quad, slice(0, 50, 2), 1)

im = ax.imshow(B_ratio, extent=[0,5,-5,5],  norm=colors.LogNorm(), cmap = 'plasma')
qr = ax.quiver(x_quad, z_quad, u, w, pivot = 'mid', scale = 30)
#fig.colorbar(im)
#plt.show()

num_range = np.linspace(0, 2*np.pi, 10)

def animate_im(num, x1, z1, a, args, R):
    phi = num
    x_quad, y_quad, z_quad, u_quad, v_quad, w_quad = B_aligned_cart(x1, z1, a, args, 2, R, phi)
    x_dip, y_dip, z_dip, u_dip, v_dip, w_dip = B_aligned_cart(x1, z1, a, args, 1, R, phi)
    
    uq_d = u_quad - u_dip
    vq_d = v_quad - v_dip
    wq_d = w_quad - w_dip
    top = np.sqrt((uq_d * uq_d) + (vq_d * vq_d) + (wq_d * wq_d))
    bottom = np.sqrt((u_dip * u_dip) + (v_dip * v_dip) + (w_dip * w_dip))
    
    B_ratio = abs(top) / abs(bottom)
    
    mag_vec = np.sqrt((u_quad * u_quad) + (w_quad * w_quad))
    u = u_quad / mag_vec
    w = w_quad / mag_vec
    
    u = np.delete(u, slice(0, 100, 2), 0)
    u = np.delete(u, slice(0, 100, 2), 1)
    u = np.delete(u, slice(0, 50, 2), 0)
    u = np.delete(u, slice(0, 50, 2), 1)

    w = np.delete(w, slice(0, 100, 2), 0)
    w = np.delete(w, slice(0, 100, 2), 1)
    w = np.delete(w, slice(0, 50, 2), 0)
    w = np.delete(w, slice(0, 50, 2), 1)
    
    im.set_array(B_ratio)
    qr.set_UVC(u, w)
    plt.title('Longitude = {:.0f}'.format(num * 180 / np.pi))
    plt.xlabel('x/a')
    plt.ylabel('z/a')
    return im, qr

anim3 = animation.FuncAnimation(fig, animate_im, frames = num_range, \
                                fargs=(x1, z1, a, args, R),interval=500, blit=False)

fig.colorbar(im)
plt.show()

#anim3.save('Harry/ratio_animation_with_vectors.gif')

#%%

z1 = np.linspace(-5 * a, 5 * a, 10)
x1 = np.linspace(0, 5 * a, 10)

phi = np.linspace(0, 2 * np.pi, 10)

fig, axs = plt.subplots(2, 5)

for i in range(len(phi)):
    
    L = np.array([2, 3, 4])
    for j in L:
        
        theta_thr = np.arcsin(j ** (-0.5)) 
    
        theta_in = np.linspace(theta_thr, np.pi - theta_thr, 100)
        
        r_L = j * np.sin(theta_in) * np.sin(theta_in)
        
        y = r_L * np.cos(theta_in)
        x = r_L * np.sin(theta_in)
        axs[int(i / 5)][i%5].plot(x, y, label = 'L = {}'.format(j))
    
    axs[int(i / 5)][i%5].legend()
    
    x_quad, y_quad, z_quad, u_quad, v_quad, w_quad = B_aligned_cart(x1, z1, a, args, 2, R, phi[i])
    x_dip, y_dip, z_dip, u_dip, v_dip, w_dip = B_aligned_cart(x1, z1, a, args, 1, R, phi[i])
    
    uq_d = u_quad - u_dip
    vq_d = v_quad - v_dip
    wq_d = w_quad - w_dip
    top = np.sqrt((uq_d * uq_d) + (vq_d * vq_d) + (wq_d * wq_d))
    bottom = np.sqrt((u_dip * u_dip) + (v_dip * v_dip) + (w_dip * w_dip))
    
    B_ratio = abs(top) / abs(bottom)
    
    mag_vec = np.sqrt((u_quad * u_quad) + (w_quad * w_quad))
    u = u_quad / mag_vec
    w = w_quad / mag_vec
    
    u = np.delete(u, slice(0, 100, 2), 0)
    u = np.delete(u, slice(0, 100, 2), 1)
    u = np.delete(u, slice(0, 50, 2), 0)
    u = np.delete(u, slice(0, 50, 2), 1)
    
    w = np.delete(w, slice(0, 100, 2), 0)
    w = np.delete(w, slice(0, 100, 2), 1)
    w = np.delete(w, slice(0, 50, 2), 0)
    w = np.delete(w, slice(0, 50, 2), 1)
    
    x_quad = np.delete(x_quad, slice(0, 100, 2), 0)
    x_quad = np.delete(x_quad, slice(0, 100, 2), 1)
    x_quad = np.delete(x_quad, slice(0, 50, 2), 0)
    x_quad = np.delete(x_quad, slice(0, 50, 2), 1)
    
    z_quad = np.delete(z_quad, slice(0, 100, 2), 0)
    z_quad = np.delete(z_quad, slice(0, 100, 2), 1)
    z_quad = np.delete(z_quad, slice(0, 50, 2), 0)
    z_quad = np.delete(z_quad, slice(0, 50, 2), 1)
    
    im = axs[int(i / 5)][i%5].imshow(B_ratio, extent=[0,5,-5,5], norm=colors.LogNorm(vmin=0.1, vmax=1.3), cmap = 'plasma')
    qr = axs[int(i / 5)][i%5].quiver(x_quad, z_quad, u, w, pivot = 'mid', scale = 30)
    axs[int(i / 5)][i%5].set_title('Longitude = {:.0f} deg'.format(phi[i] * 180 / np.pi))
    axs[int(i / 5)][i%5].set(xlabel='x / a', ylabel='z / a')
    fig.colorbar(im, ax= axs[int(i / 5)][i%5])
    
plt.show()
#%%

z1 = np.linspace(-5 * a, 5 * a, 20)
x1 = np.linspace(0, 5 * a, 10)

phi = np.linspace(0, 2 * np.pi, 100)

B_all = []

for i in tqdm(phi):

    x_quad, y_quad, z_quad, u_quad, v_quad, w_quad = B_aligned_cart(x1, z1, a, args, 2, R, i)
    x_dip, y_dip, z_dip, u_dip, v_dip, w_dip = B_aligned_cart(x1, z1, a, args, 1, R, i)
    
    uq_d = u_quad - u_dip
    vq_d = v_quad - v_dip
    wq_d = w_quad - w_dip
    top = np.sqrt((uq_d * uq_d) + (vq_d * vq_d) + (wq_d * wq_d))
    bottom = np.sqrt((u_dip * u_dip) + (v_dip * v_dip) + (w_dip * w_dip))
    
    # B_mag_quad = np.nan_to_num(B_mag_quad)
    # B_mag_dip = np.nan_to_num(B_mag_dip)
    
    B_ratio = abs(top) / abs(bottom)
        
    B_all.append(np.nanmax(B_ratio))
    
    
    
plt.plot(phi * 180 / np.pi, B_all)
plt.xlabel('Longitude (deg)')
plt.ylabel('Max ratio')
plt.show()

#%%
L = 2

theta_thr = np.arcsin(L ** (-0.5)) 

theta_in = np.linspace(theta_thr, np.pi - theta_thr, 100)
phi_in = np.linspace(0, 2 * np.pi, 100)
theta_deg = theta_in * 180 / np.pi
phi_deg = phi_in * 180 / np.pi

x_quad, y_quad, z_quad, u_quad, v_quad, w_quad = B_Lshell(L, theta_in, phi_in, a, args, 2, R)
x_dip, y_dip, z_dip, u_dip, v_dip, w_dip = B_Lshell(L, theta_in, phi_in, a, args, 1, R)

uq_d = u_quad - u_dip
vq_d = v_quad - v_dip
wq_d = w_quad - w_dip
top = np.sqrt((uq_d * uq_d) + (vq_d * vq_d) + (wq_d * wq_d))
bottom = np.sqrt((u_dip * u_dip) + (v_dip * v_dip) + (w_dip * w_dip))

# B_mag_quad = np.nan_to_num(B_mag_quad)
# B_mag_dip = np.nan_to_num(B_mag_dip)

B_ratio = abs(top) / abs(bottom)

B_ratio = B_ratio.reshape((len(theta_in), len(phi_in)))

u_all = []
v_all = []

for i in range(len(u_quad)):
    uvw_quad = np.array([u_quad[i], v_quad[i], w_quad[i]])
    quad_norm = uvw_quad / np.linalg.norm(uvw_quad)
    
    uvw_dip = np.array([u_dip[i], v_dip[i], w_dip[i]])
    dip_norm = uvw_dip / np.linalg.norm(uvw_dip)
    
    dot = np.dot(quad_norm, dip_norm)
    angle = np.arccos(dot)

    rot_matrix = np.array([[np.cos(angle), -np.sin(angle)],
                       [np.sin(angle), np.cos(angle)]])
    
    u, v = np.matmul(rot_matrix, np.array([0, 1]))
    u_all.append(u)
    v_all.append(v)
    
u = np.array(u_all)
v = np.array(v_all)
u = u.reshape((len(theta_in), len(phi_in)))    
v = v.reshape((len(theta_in), len(phi_in)))

u = np.delete(u, slice(0, 100, 2), 0)
u = np.delete(u, slice(0, 100, 2), 1)
u = np.delete(u, slice(0, 50, 2), 0)
u = np.delete(u, slice(0, 50, 2), 1)

v = np.delete(v, slice(0, 100, 2), 0)
v = np.delete(v, slice(0, 100, 2), 1)
v = np.delete(v, slice(0, 50, 2), 0)
v = np.delete(v, slice(0, 50, 2), 1)

phi_vec = phi_deg[::2]
phi_vec = phi_vec[::2]

theta_vec = theta_deg[::2]
theta_vec = theta_vec[::2]

fig, ax = plt.subplots(1, 1)
im = ax.imshow(B_ratio,  extent=[phi_deg[0],phi_deg[-1],theta_deg[-1],theta_deg[0]], cmap = 'plasma', norm=colors.LogNorm(), aspect='auto')
fig.colorbar(im)
ax.quiver(phi_vec, theta_vec, u, v, scale = 30, pivot = 'mid')
plt.xlabel('Longitude (deg)', fontsize=16)
plt.ylabel('Latitude (deg)', fontsize=16)
plt.title('Max ratio, L = 2', fontsize=16)
plt.show()
#%%
L = np.arange(2, 10)
fig, axs = plt.subplots(2, 4)

for j in range(len(L)):

    theta_thr = np.arcsin(L[j] ** (-0.5)) 
    
    theta_in = np.linspace(theta_thr, np.pi - theta_thr, 100)
    phi_in = np.linspace(0, 2 * np.pi, 100)
    theta_deg = theta_in * 180 / np.pi
    phi_deg = phi_in * 180 / np.pi
    
    x_quad, y_quad, z_quad, u_quad, v_quad, w_quad = B_Lshell(L[j], theta_in, phi_in, a, args, 2, R)
    x_dip, y_dip, z_dip, u_dip, v_dip, w_dip = B_Lshell(L[j], theta_in, phi_in, a, args, 1, R)
    
    uq_d = u_quad - u_dip
    vq_d = v_quad - v_dip
    wq_d = w_quad - w_dip
    top = np.sqrt((uq_d * uq_d) + (vq_d * vq_d) + (wq_d * wq_d))
    bottom = np.sqrt((u_dip * u_dip) + (v_dip * v_dip) + (w_dip * w_dip))
    
    B_ratio = abs(top) / abs(bottom)
    
    B_ratio = B_ratio.reshape((len(theta_in), len(phi_in)))
    
    u_all = []
    v_all = []
    
    for i in range(len(u_quad)):
        uvw_quad = np.array([u_quad[i], v_quad[i], w_quad[i]])
        quad_norm = uvw_quad / np.linalg.norm(uvw_quad)
        
        uvw_dip = np.array([u_dip[i], v_dip[i], w_dip[i]])
        dip_norm = uvw_dip / np.linalg.norm(uvw_dip)
        
        dot = np.dot(quad_norm, dip_norm)
        angle = np.arccos(dot)
    
        rot_matrix = np.array([[np.cos(angle), -np.sin(angle)],
                           [np.sin(angle), np.cos(angle)]])
        
        u, v = np.matmul(rot_matrix, np.array([0, 1]))
        u_all.append(u)
        v_all.append(v)
        
    u = np.array(u_all)
    v = np.array(v_all)
    u = u.reshape((len(theta_in), len(phi_in)))    
    v = v.reshape((len(theta_in), len(phi_in)))
    
    u = np.delete(u, slice(0, 100, 2), 0)
    u = np.delete(u, slice(0, 100, 2), 1)
    u = np.delete(u, slice(0, 50, 2), 0)
    u = np.delete(u, slice(0, 50, 2), 1)
    
    v = np.delete(v, slice(0, 100, 2), 0)
    v = np.delete(v, slice(0, 100, 2), 1)
    v = np.delete(v, slice(0, 50, 2), 0)
    v = np.delete(v, slice(0, 50, 2), 1)
    
    phi_vec = phi_deg[::2]
    phi_vec = phi_vec[::2]
    
    theta_vec = theta_deg[::2]
    theta_vec = theta_vec[::2]
    
    im = axs[int(j / 4)][j%4].imshow(B_ratio,  extent=[phi_deg[0],phi_deg[-1],\
            theta_deg[-1],theta_deg[0]], cmap = 'plasma', norm=colors.LogNorm(vmin = 0.1, vmax = 1.369632402682142), aspect='auto')
    axs[int(j / 4)][j%4].quiver(phi_vec, theta_vec, u, v, scale = 40, headwidth = 6, pivot = 'mid')
    axs[int(j / 4)][j%4].set_title('Max ratio, L = {}'.format(L[j]))
    axs[int(j / 4)][j%4].set(xlabel='Longitude (deg)', ylabel='Latitude (deg)')
    fig.colorbar(im, ax= axs[int(j / 4)][j%4])
plt.show()
#%%
# z1 = np.linspace(-5 * a, 5 * a, 20)
# x1 = np.linspace(0, 5 * a, 10)
# y1 = [0]

# coords, dB, B = gradB_fun(x1, y1, z1, a, g, h, 2, R)

#test points

x1 = [0 * a]
y1 = [0 * a]
z1 = [4 * a]

coords, dB1, B = gradB_fun(x1, y1, z1, a, g, h, 2, R)

print('A',dB1)

x1 = [1 * a]
y1 = [2 * a]
z1 = [3 * a]

coords, dB2, B = gradB_fun(x1, y1, z1, a, g, h, 2, R)

print('B', dB2, np.linalg.norm(dB2))

x1 = [-4 * a]
y1 = [2.5 * a]
z1 = [6 * a]

coords, dB3, B = gradB_fun(x1, y1, z1, a, g, h, 2, R)

print('C', dB3)

#%%
x1 = np.linspace(0, 5 * a, 10)
y1 = [0]
z1 = np.linspace(-5 * a, 5 * a, 20)

coords, dB, B = gradB_fun(x1, y1, z1, a, g, h, 1, R)

du = dB[0]
dv = dB[1]
dw = dB[2]

du = du.reshape((len(x1), len(z1))).T
dv = dv.reshape((len(x1), len(z1))).T
dw = dw.reshape((len(x1), len(z1))).T

x, y, z = coords

x = x.reshape((len(x1), len(z1))).T
y = y.reshape((len(x1), len(z1))).T
z = z.reshape((len(x1), len(z1))).T

plt.quiver(x, z, du, dw)

#%%

x1 = np.linspace(0, 5 * a, 20)
z1 = np.linspace(-5 * a, 5 * a, 20)
phi_rot = 1 * np.pi#np.linspace(0, 2 * np.pi, 5)

coords, dB, B = gradB_long(x1, z1, a, g, h, 1, R, phi_rot)

x, y, z = coords
du, dv, dw = dB
u, v, w = B

plt.quiver(x, z, du, dw, pivot = 'mid')
plt.show()

#%%
x1 = np.linspace(0, 5 * a, 100)
z1 = np.linspace(-5 * a, 5 * a, 100)
phi_rot = 1.5 * np.pi#np.linspace(0, 2 * np.pi, 5)

coords, dB, B = gradB_long(x1, z1, a, g, h, 2, R, phi_rot)

x, y, z = coords
du, dv, dw = dB
u, v, w = B

drift = np.zeros((3, len(u), len(u[0])))

for i in range(len(u)):
    for j in range(len(u[0])):
        dB_temp = np.array([du[i][j], dv[i][j], dw[i][j]])
        B_temp = np.array([u[i][j], v[i][j], w[i][j]])
        
        cross_temp = np.cross(B_temp, dB_temp)
        drift[0][i][j] = cross_temp[0] / np.linalg.norm(cross_temp)
        drift[1][i][j] = cross_temp[1] / np.linalg.norm(cross_temp)
        drift[2][i][j] = cross_temp[2] / np.linalg.norm(cross_temp)

u_dr, v_dr, w_dr = drift

u_dr = np.delete(u_dr, slice(0, 100, 2), 0)
u_dr = np.delete(u_dr, slice(0, 100, 2), 1)
u_dr = np.delete(u_dr, slice(0, 50, 2), 0)
u_dr = np.delete(u_dr, slice(0, 50, 2), 1)

w_dr = np.delete(w_dr, slice(0, 100, 2), 0)
w_dr = np.delete(w_dr, slice(0, 100, 2), 1)
w_dr = np.delete(w_dr, slice(0, 50, 2), 0)
w_dr = np.delete(w_dr, slice(0, 50, 2), 1)

x = np.delete(x, slice(0, 100, 2), 0)
x = np.delete(x, slice(0, 100, 2), 1)
x = np.delete(x, slice(0, 50, 2), 0)
x = np.delete(x, slice(0, 50, 2), 1)

z = np.delete(z, slice(0, 100, 2), 0)
z = np.delete(z, slice(0, 100, 2), 1)
z = np.delete(z, slice(0, 50, 2), 0)
z = np.delete(z, slice(0, 50, 2), 1)

#%%
im = plt.imshow(v_dr[::-1],  extent=[0, 5, -5, 5], cmap = 'plasma', vmin = -1, vmax = 1)#, aspect='auto')
plt.quiver(x, z[::-1], u_dr, w_dr[::-1], pivot = 'mid', width = 0.005, scale = 4, scale_units = 'x')
plt.colorbar(im)
plt.xlabel('x/a', fontsize=16)
plt.ylabel('z/a', fontsize=16)
plt.title('Drift', fontsize=16)
plt.show()

#%%
x1 = np.linspace(0, 5 * a, 100)
z1 = np.linspace(-5 * a, 5 * a, 100)
phi = np.linspace(0, 2 * np.pi, 10)

fig, axs = plt.subplots(2, 5)

for k in range(len(phi)):

    coords, dB, B = gradB_long(x1, z1, a, g, h, 2, R, phi[k])
    
    x, y, z = coords
    du, dv, dw = dB
    u, v, w = B
    
    drift = np.zeros((3, len(u), len(u[0])))
    
    for i in range(len(u)):
        for j in range(len(u[0])):
            dB_temp = np.array([du[i][j], dv[i][j], dw[i][j]])
            B_temp = np.array([u[i][j], v[i][j], w[i][j]])
            
            cross_temp = np.cross(B_temp, dB_temp)
            drift[0][i][j] = cross_temp[0] / np.linalg.norm(cross_temp)
            drift[1][i][j] = cross_temp[1] / np.linalg.norm(cross_temp)
            drift[2][i][j] = cross_temp[2] / np.linalg.norm(cross_temp)
    
    u_dr, v_dr, w_dr = drift
    
    u_dr = np.delete(u_dr, slice(0, 100, 2), 0)
    u_dr = np.delete(u_dr, slice(0, 100, 2), 1)
    u_dr = np.delete(u_dr, slice(0, 50, 2), 0)
    u_dr = np.delete(u_dr, slice(0, 50, 2), 1)
    
    w_dr = np.delete(w_dr, slice(0, 100, 2), 0)
    w_dr = np.delete(w_dr, slice(0, 100, 2), 1)
    w_dr = np.delete(w_dr, slice(0, 50, 2), 0)
    w_dr = np.delete(w_dr, slice(0, 50, 2), 1)
    
    x = np.delete(x, slice(0, 100, 2), 0)
    x = np.delete(x, slice(0, 100, 2), 1)
    x = np.delete(x, slice(0, 50, 2), 0)
    x = np.delete(x, slice(0, 50, 2), 1)
    
    z = np.delete(z, slice(0, 100, 2), 0)
    z = np.delete(z, slice(0, 100, 2), 1)
    z = np.delete(z, slice(0, 50, 2), 0)
    z = np.delete(z, slice(0, 50, 2), 1)
    
    axs[int(k / 5)][k%5].quiver(x, z[::-1], u_dr, w_dr[::-1], pivot = 'mid', width = 0.005, scale = 4, scale_units = 'x')
    im = axs[int(k / 5)][k%5].imshow(v_dr[::-1],  extent=[0, 5, -5, 5], cmap = 'plasma', vmin = -1, vmax = 1)#, aspect='auto')
    fig.colorbar(im, ax= axs[int(k / 5)][k%5])
    axs[int(k / 5)][k%5].set(xlabel='x/a', ylabel='z/a')
    axs[int(k / 5)][k%5].set_title('Longitude = {:.0f} deg'.format(phi[k] * 180 / np.pi))
plt.show()

#%%
L = 3

theta_thr = np.arcsin(L ** (-0.5)) 

theta_in = np.linspace(theta_thr, np.pi - theta_thr, 15)
phi_in = [np.pi * 1]#np.linspace(0, 2 * np.pi, 10)

for i in theta_in:
    # theta_deg = theta_in * 180 / np.pi
    # phi_deg = phi_in * 180 / np.pi
    
    coords, dB, B, normal, phi_hat, sigma = dB_Lshell(L, [i], phi_in, a, g, h, 2, R)
    
    du, dv, dw = dB
    x, y, z = coords
    x = x[0]
    y = y[0]
    z = z[0]
    
    
    phi_rot = phi_in[0]
    
    R_z = np.array([[np.cos(phi_rot), - np.sin(phi_rot), 0],
                   [np.sin(phi_rot), np.cos(phi_rot), 0],
                   [0, 0, 1]])
    
    normal = np.matmul(np.linalg.inv(R_z), normal[0])
    coords = np.matmul(np.linalg.inv(R_z), coords)

    phi_hat = np.matmul(np.linalg.inv(R_z), phi_hat[0])
    
    sigma = np.matmul(np.linalg.inv(R_z), sigma[0])
    
    t = np.linspace(0, np.pi, 1000)
    p = np.array([0]) 
    r = 3 * np.sin(t) * np.sin(t)
    
    xp, yp, zp = Sph_to_Cart(r, t, p)
    
    plt.quiver(coords[0], coords[2], normal[0], normal[2])
    plt.quiver(coords[0], coords[2], sigma[0], sigma[2])
plt.plot(xp, zp)
plt.ylim(-1.5, 1.5)
plt.xlabel('x/a', fontsize = 16)
plt.ylabel('z/a', fontsize = 16)
plt.title('Direction test')
plt.show()

#%%
L = 3

theta_thr = np.arcsin(L ** (-0.5)) 

theta_in = np.linspace(theta_thr, np.pi - theta_thr, 100)
phi_in = np.linspace(0, 2 * np.pi, 100)
theta_deg = theta_in * 180 / np.pi
phi_deg = phi_in * 180 / np.pi

coords, dB, B, normal, phi_hat, sigma = dB_Lshell(L, theta_in, phi_in, a, g, h, 2, R)

du, dv, dw = dB
x, y, z = coords

drift_n = []
drift_sigma = []
drift_phi = []
# b_n = []
# b_sigma = []
# b_phi = []
phi_vec = []
theta_vec = []

for i in range(len(dB[0])):
    Bt = np.array([B[0][i], B[1][i], B[2][i]])
    dBt = np.array([dB[0][i], dB[1][i], dB[2][i]])
    dt = np.cross(Bt, dBt)
    dt = dt / np.linalg.norm(dt)
    
    drift_n.append(np.dot(dt, normal[i]))
    drift_sigma.append(np.dot(dt, sigma[i]))
    drift_phi.append(np.dot(dt, phi_hat[i]))
    
    # b_n.append(np.dot(Bt, normal[i]))
    # b_sigma.append(np.dot(Bt, sigma[i]))
    # b_phi.append(np.dot(Bt, phi_hat[i]))

drift_n = np.array(drift_n).reshape((len(theta_in), len(phi_in)))
drift_sigma = np.array(drift_sigma).reshape((len(theta_in), len(phi_in)))
drift_phi = np.array(drift_phi).reshape((len(theta_in), len(phi_in)))

# b_n = np.array(b_n).reshape((len(theta_in), len(phi_in)))
# b_sigma = np.array(b_sigma).reshape((len(theta_in), len(phi_in)))
# b_phi = np.array(b_phi).reshape((len(theta_in), len(phi_in)))

drift_phi = np.delete(drift_phi, slice(0, 100, 2), 0)
drift_phi = np.delete(drift_phi, slice(0, 100, 2), 1)
drift_phi = np.delete(drift_phi, slice(0, 50, 2), 0)
drift_phi = np.delete(drift_phi, slice(0, 50, 2), 1)

drift_sigma = np.delete(drift_sigma, slice(0, 100, 2), 0)
drift_sigma = np.delete(drift_sigma, slice(0, 100, 2), 1)
drift_sigma = np.delete(drift_sigma, slice(0, 50, 2), 0)
drift_sigma = np.delete(drift_sigma, slice(0, 50, 2), 1)

phi_vec = phi_deg[::2]
phi_vec = phi_vec[::2]

theta_vec = theta_deg[::2]
theta_vec = theta_vec[::2]

im = plt.imshow(drift_n,  extent=[phi_deg[0],phi_deg[-1],\
        theta_deg[-1],theta_deg[0]], cmap = 'plasma', aspect = 'auto')
plt.quiver(phi_vec, theta_vec, drift_phi, drift_sigma, pivot = 'mid')
plt.colorbar(im)
#.quiver(phi_deg, theta_deg, b_phi, b_sigma, pivot = 'mid')
plt.xlabel('Longitude (deg)', fontsize = 16)
plt.ylabel('Latitude (deg)', fontsize = 16)
plt.title('L = 3', fontsize = 16)
plt.show()

#%%
L = np.arange(2, 10)
fig, axs = plt.subplots(2, 4)

for j in range(len(L)):

    theta_thr = np.arcsin(L[j] ** (-0.5))
    
    theta_in = np.linspace(theta_thr, np.pi - theta_thr, 100)
    phi_in = np.linspace(0, 2 * np.pi, 100)
    theta_deg = theta_in * 180 / np.pi
    phi_deg = phi_in * 180 / np.pi

    coords, dB, B, normal, phi_hat, sigma = dB_Lshell(L[j], theta_in, phi_in, a, g, h, 2, R)

    du, dv, dw = dB
    x, y, z = coords

    drift_n = []
    drift_sigma = []
    drift_phi = []
    phi_vec = []
    theta_vec = []

    for i in range(len(dB[0])):
        Bt = np.array([B[0][i], B[1][i], B[2][i]])
        dBt = np.array([dB[0][i], dB[1][i], dB[2][i]])
        dt = np.cross(Bt, dBt)
        dt = dt / np.linalg.norm(dt)
        
        drift_n.append(np.dot(dt, normal[i]))
        drift_sigma.append(np.dot(dt, sigma[i]))
        drift_phi.append(np.dot(dt, phi_hat[i]))

    drift_n = np.array(drift_n).reshape((len(theta_in), len(phi_in)))
    drift_sigma = np.array(drift_sigma).reshape((len(theta_in), len(phi_in)))
    drift_phi = np.array(drift_phi).reshape((len(theta_in), len(phi_in)))

    drift_phi = np.delete(drift_phi, slice(0, 100, 2), 0)
    drift_phi = np.delete(drift_phi, slice(0, 100, 2), 1)
    drift_phi = np.delete(drift_phi, slice(0, 50, 2), 0)
    drift_phi = np.delete(drift_phi, slice(0, 50, 2), 1)

    drift_sigma = np.delete(drift_sigma, slice(0, 100, 2), 0)
    drift_sigma = np.delete(drift_sigma, slice(0, 100, 2), 1)
    drift_sigma = np.delete(drift_sigma, slice(0, 50, 2), 0)
    drift_sigma = np.delete(drift_sigma, slice(0, 50, 2), 1)

    phi_vec = phi_deg[::2]
    phi_vec = phi_vec[::2]

    theta_vec = theta_deg[::2]
    theta_vec = theta_vec[::2]

    im = axs[int(j / 4)][j%4].imshow(drift_n,  extent=[phi_deg[0],phi_deg[-1],\
            theta_deg[-1],theta_deg[0]], cmap = 'plasma', aspect='auto')
    axs[int(j / 4)][j%4].quiver(phi_vec, theta_vec, drift_phi, drift_sigma, scale = 40, headwidth = 6, pivot = 'mid')
    axs[int(j / 4)][j%4].set_title('Max ratio, L = {}'.format(L[j]))
    axs[int(j / 4)][j%4].set(xlabel='Longitude (deg)', ylabel='Latitude (deg)')
    fig.colorbar(im, ax= axs[int(j / 4)][j%4])
    #fig.suptitle('Dipole only', fontsize=16)
plt.show()
#%%

#loadBField('Output/complete_field_phi=0_CI.npz')
#%%

xn, yn, zn, un, vn, wn = loadBField('Output/dipole_nick_fieldaligned.npz')

xc, yc, zc, uc, vc, wc = reshapeCN(x, y, z, u, v, w, 4, 4, 4)

#%%

print(un - uc)

#print(x, y, z, u, v, w)
