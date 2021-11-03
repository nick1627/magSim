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

args = np.array([g01, g11, h11, g02, g12, h12, g22, h22])

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
x1 = np.linspace(-1.5 * a, 1.5 * a, 4)
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
plt.rcParams["figure.figsize"] = [10, 8]
plt.rcParams["figure.autolayout"] = True
plt.scatter(x=x,y=z,c=B_ratio, s=40, norm=matplotlib.colors.LogNorm())
plt.colorbar()
#plt.contourf(X_mat, Z_mat, B_mag_mat)
plt.ylabel('$z / \u03B1$')
plt.quiver(x[::13], z[::13], u[::13], w[::13], scale = 1000000, pivot = 'mid')
plt.show()

#%%

r = np.linspace(1*a, 2 * a, 10)
theta = np.linspace(0, np.pi, 10)
phi = np.linspace(0, 2*np.pi, 20)#np.array([0, np.pi])

B_rall = Get_maxB_ratio(r, theta, phi, a, args, R, True) 
    
#     plt.scatter(x=x,y=z,c=B_ratio, s=150)#, norm=matplotlib.colors.LogNorm())
#     plt.colorbar()
#     plt.quiver(x[::13], z[::13], u[::13], w[::13], scale = 1000000, pivot = 'mid')
#     plt.title('Longitude = {:.1f}'.format(i  * 180 / np.pi))
    
#     plt.show()

#%%

r = np.linspace(1*a, 2 * a, 100)
theta = np.linspace(0, np.pi, 100)
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
z1 = np.linspace(-5 * a, 5 * a, 100)
x1 = np.linspace(0, 5 * a, 100)

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

im = ax.imshow(B_ratio, extent=[0,5,-5,5],  norm=colors.LogNorm())
fig.colorbar(im)
plt.show()
#plt.quiver(x, z, u, w)
#ax.colorbar()

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

    im.set_array(B_ratio)
    plt.title('Longitude = {:.0f}'.format(num * 180 / np.pi))
    plt.xlabel('x/a')
    plt.ylabel('z/a')
    return [im]

anim3 = animation.FuncAnimation(fig, animate_im, frames = num_range, \
                                fargs=(x1, z1, a, args, R),interval=500, blit=False)

fig.colorbar(im)
plt.show()
#%%
#anim3.save('Harry/ratio_animation.gif')
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
        
    B_all.append(max(B_ratio))

plt.plot(phi * 180 / np.pi, B_all)
plt.xlabel('Longitude (deg)')
plt.ylabel('Max ratio')
plt.show()
#%%
#plt.scatter(x, z, c = B_mag, s = 150)
phi = 137 * np.pi / 180

x,y,z,u,v,w = B_mag_cart(x1, z1, a, args, 2, R, phi,True)
saveBField(x, y, z, u, v, w, 'Output/complete_field_phi=137_CI')
#%%

#loadBField('Output/complete_field_phi=0_CI.npz')
#%%

xn, yn, zn, un, vn, wn = loadBField('Output/dipole_nick_fieldaligned.npz')

xc, yc, zc, uc, vc, wc = reshapeCN(x, y, z, u, v, w, 4, 4, 4)

#%%

print(un - uc)

#print(x, y, z, u, v, w)
