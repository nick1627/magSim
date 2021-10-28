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

#%%
#CONSTANTS    
a = 25600000 # Uranus' radius

#%%
#Uranus' coefficients
g01 = 1#11278 #z
g11 = 1#10928 #x
h11 = 1#-16049 #y
g02 = 1#-9648
g12 = 1#-12284
h12 = 1#6405
g22 = 1#1453
h22 = 1#4220

args = np.array([g01, g11, h11, g02, g12, h12, g22, h22])

r = np.linspace(2 * a, 3 * a, 5)
theta = np.linspace(0, np.pi , 10)
phi = np.linspace(0, 2 * np.pi, 10)

x, y, z, u, v, w = Get_B_sph(r, theta, phi, a, args, 2)

ThreeD_plot(x, y, z, u, v, w)

#%%
z1 = [0]#np.linspace(-1.5 * a, 1.5 * a, 4)
y1 = np.linspace(-1.5 * a, 1.5 * a, 20)
x1 = np.linspace(-1.5 * a, 1.5 * a, 20)

g01 = 11278 #z
g11 = 10928 #x
h11 = -16049 #y
g02 = -9648
g12 = -12284
h12 = 6405
g22 = 1453
h22 = 4220

args = np.array([g01, g11, h11, g02, g12, h12, g22, h22])

x, y, z, u, v, w = Get_B_cart(x1, y1, z1, a, args, 2, True)

TwoD_plot(x, y, u, v, 'z')
            
#%%
x1 = np.linspace(-1.5 * a, 1.5 * a, 4)
y1 = np.linspace(-1.5 * a, 1.5 * a, 4)
z1 = np.linspace(-1.5 * a, 1.5 * a, 4)
    
    
g01 = 1#11278 #z
g11 = 0#10928 #x
h11 = 0#-16049 #y
g02 = 0#-9648
g12 = 0#-12284
h12 = 0#6405
g22 = 0#1453
h22 = 0#4220


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

g01 = 22454.1 #z
g11 = 0#10928 #x
h11 = 0#-16049 #y
g02 = 0#-9648
g12 = 0#-12284
h12 = 0#6405
g22 = 0#1453
h22 = 0#4220

args = np.array([g01, g11, h11, g02, g12, h12, g22, h22])

x, y, z, u, v, w = Get_B_cart(x1, y1, z1, a, args, 2, True)

TwoD_plot(x, y, u, v, 'z')

#%%
x1 = np.linspace(-2* a, 2 * a, 5)
z1 = np.linspace(-2 * a, 2 * a, 5)
y1 = np.linspace(-2 * a, 2 * a, 5)

g01 = 11278 #z
g11 = 10928 #x
h11 = -16049 #y
g02 = 0#-9648
g12 = 0#-12284
h12 = 0#6405
g22 = 0#1453
h22 = 0#4220

args = np.array([g01, g11, h11, g02, g12, h12, g22, h22])

x, y, z, u, v, w = Get_B_cart(x1, y1, z1, a, args, 2)

ThreeD_plot(x, y, z, u, v, w)

#%%
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
#%%

g01 = 11278 #z
g11 = 10928 #x
h11 = -16049 #y
g02 = -9648
g12 = -12284
h12 = 6405
g22 = 1453
h22 = 4220

# y1 = [0]#np.linspace(-1.5 * a, 1.5 * a, 4)
# x1 = np.linspace(-1.5 * a, 1.5 * a, 20)
# z1 = np.linspace(-1.5 * a, 1.5 * a, 20)

plt.rcParams["figure.figsize"] = [8, 8]
plt.rcParams["figure.autolayout"] = True

r = np.linspace(a, 2 * a, 4)
theta = np.linspace(0, np.pi, 10)
phi = np.array([0, np.pi])

args = np.array([g01, g11, h11, g02, g12, h12, g22, h22])

# x, y, z, u, v, w = Get_B_cart_rot(x1, y1, z1, a, args, 2, R)
x, y, z, u, v, w = Get_B_sph_rot(r, theta, phi, a, args, 2, R)

R_z = np.array([[np.cos(-phi[0]), - np.sin(-phi[0]), 0],
               [np.sin(-phi[0]), np.cos(-phi[0]), 0],
               [0, 0, 1]])

xyz = np.array([x, y, z])
x, y, z = np.matmul(R_z, xyz)
uvw = np.array([u, v, w])
u, v, w = np.matmul(R_z, uvw)

fig, ax = plt.subplots(1, 1)
qr = ax.quiver(x, z, u, w, color='black', pivot = 'mid')
plt.title('Phi = 0')

num_range = np.linspace(0, np.pi, 5)

def animate(num, qr, r, theta):
   phi = np.array([num, num + np.pi])

   args = np.array([g01, g11, h11, g02, g12, h12, g22, h22])

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

anim = animation.FuncAnimation(fig, animate, frames = num_range, fargs=(qr, r, theta),
                              interval=1000, blit=False)

Circle1 = plt.Circle((0, 0), 1, color = 'red', fill = False)
ax.add_patch(Circle1)

plt.show()
#TwoD_plot(x, z, u, w, 'y')

anim.save('Harry/quadrupole_test.gif')

#%%

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
# y1 = [0]#np.linspace(-1.5 * a, 1.5 * a, 4)
# x1 = np.linspace(-1.5 * a, 1.5 * a, 20)
# z1 = np.linspace(-1.5 * a, 1.5 * a, 20)

r = np.linspace(a, 2 * a, 100)
theta = np.linspace(0, np.pi, 100)
phi = np.linspace(0, 2*np.pi, 100)#np.array([0, np.pi])

# x, y, z, u, v, w = Get_B_cart_rot(x1, y1, z1, a, args, 2, R)
x, y, z, u, v, w = Get_B_sph_rot(r, theta, phi, a, args, 2, R)
x2, y2, z2, u2, v2, w2 = Get_B_sph_rot(r, theta, phi, a, args, 1, R)

R_z = np.array([[np.cos(-phi[0]), - np.sin(-phi[0]), 0],
               [np.sin(-phi[0]), np.cos(-phi[0]), 0],
               [0, 0, 1]])

xyz = np.array([x, y, z])
x, y, z = np.matmul(R_z, xyz)
uvw = np.array([u, v, w])
u, v, w = np.matmul(R_z, uvw)

xyz2 = np.array([x2, y2, z2])
x2, y2, z2 = np.matmul(R_z, xyz2)
uvw2 = np.array([u2, v2, w2])
u2, v2, w2 = np.matmul(R_z, uvw2)

B_mag = np.sqrt((u * u) + (v * v) + (w * w))
B_mag2 = np.sqrt((u2 * u2) + (v2 * v2) + (w2 * w2))

B_ratio = B_mag2 / B_mag
# # X_mat = np.zeros((8, 10))
# # Z_mat = np.zeros((8, 10))
# # B_mag_mat = np.zeros((8,10))


# # for i in range(len(x)):
#     s = int(i/10)
#     d = i%10
#     X_mat[s][d] = x[i]
#     Z_mat[s][d] = z[i]
#     B_mag_mat[s][d] = B_mag[i]
# B_mag = B_mag / max(B_mag)
B_ratio = B_ratio / max(B_ratio)

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

r = np.linspace(1*a, 2 * a, 100)
theta = np.linspace(0, np.pi, 100)
phi = np.linspace(0, 2*np.pi, 1)#np.array([0, np.pi])

args = np.array([g01, g11, h11, g02, g12, h12, g22, h22])

B_all = []

for i in tqdm(phi):

    # x, y, z, u, v, w = Get_B_cart_rot(x1, y1, z1, a, args, 2, R)
    x, y, z, u, v, w = Get_B_sph_rot(r, theta, [i], a, args, 2, R)
    x2, y2, z2, u2, v2, w2 = Get_B_sph_rot(r, theta, [i], a, args, 1, R)
    
    R_z = np.array([[np.cos(-i), - np.sin(-i), 0],
                   [np.sin(-i), np.cos(-i), 0],
                   [0, 0, 1]])
    
    xyz = np.array([x, y, z])
    x, y, z = np.matmul(R_z, xyz)
    uvw = np.array([u, v, w])
    u, v, w = np.matmul(R_z, uvw)
    
    xyz2 = np.array([x2, y2, z2])
    x2, y2, z2 = np.matmul(R_z, xyz2)
    uvw2 = np.array([u2, v2, w2])
    u2, v2, w2 = np.matmul(R_z, uvw2)
    
    B_mag = np.sqrt((u * u) + (v * v) + (w * w))
    B_mag2 = np.sqrt((u2 * u2) + (v2 * v2) + (w2 * w2))
    
    B_ratio = abs((B_mag - B_mag2) / B_mag2)
    
    B_all.append(max(B_ratio))
    
    plt.scatter(x=x,y=z,c=B_ratio, s=150)#, norm=matplotlib.colors.LogNorm())
    plt.colorbar()
    plt.quiver(x[::13], z[::13], u[::13], w[::13], scale = 1000000, pivot = 'mid')
    plt.title('Longitude = {:.1f}'.format(i  * 180 / np.pi))
    
    plt.show()
#%%

plt.plot(phi * 180 / np.pi, B_all)
plt.xlabel('Longitude')
plt.ylabel('Max ratio')

#%%

r = np.linspace(1*a, 2 * a, 100)
theta = np.linspace(0, np.pi, 100)
phi = [0]#np.linspace(0, 2*np.pi, 50)#np.array([0, np.pi])

x, y, z, u, v, w = Get_B_sph_rot(r, theta, [0], a, args, 2, R)
x2, y2, z2, u2, v2, w2 = Get_B_sph_rot(r, theta, [0], a, args, 1, R)

R_z = np.array([[np.cos(0), - np.sin(0), 0],
                [np.sin(0), np.cos(0), 0],
                [0, 0, 1]])

xyz = np.array([x, y, z])
x, y, z = np.matmul(R_z, xyz)
uvw = np.array([u, v, w])
u, v, w = np.matmul(R_z, uvw)

xyz2 = np.array([x2, y2, z2])
x2, y2, z2 = np.matmul(R_z, xyz2)
uvw2 = np.array([u2, v2, w2])
u2, v2, w2 = np.matmul(R_z, uvw2)

B_mag = np.sqrt((u * u) + (v * v) + (w * w))
B_mag2 = np.sqrt((u2 * u2) + (v2 * v2) + (w2 * w2))

B_ratio = abs((B_mag - B_mag2) / B_mag2)

fig, ax = plt.subplots(1, 1)
scat = plt.scatter(x=x,y=z,c=B_ratio, s=150)#, norm=matplotlib.colors.LogNorm())
#plt.title('Phi = 0')

num_range = np.linspace(0, 2 * np.pi, 10)

def animate2(num, scat, r, theta):
   phi = [num]

   args = np.array([g01, g11, h11, g02, g12, h12, g22, h22])

   x, y, z, u, v, w = Get_B_sph_rot(r, theta, phi, a, args, 2, R)
   x2, y2, z2, u2, v2, w2 = Get_B_sph_rot(r, theta, phi, a, args, 1, R)

   R_z = np.array([[np.cos(-phi[0]), - np.sin(-phi[0]), 0],
                  [np.sin(-phi[0]), np.cos(-phi[0]), 0],
                  [0, 0, 1]])

   xyz = np.array([x, y, z])
   x, y, z = np.matmul(R_z, xyz)
   uvw = np.array([u, v, w])
   u, v, w = np.matmul(R_z, uvw)

   xyz2 = np.array([x2, y2, z2])
   x2, y2, z2 = np.matmul(R_z, xyz2)
   uvw2 = np.array([u2, v2, w2])
   u2, v2, w2 = np.matmul(R_z, uvw2)

   B_mag = np.sqrt((u * u) + (v * v) + (w * w))
   B_mag2 = np.sqrt((u2 * u2) + (v2 * v2) + (w2 * w2))

   B_ratio = abs((B_mag - B_mag2) / B_mag2)
   
   scat.set_array(B_ratio)
   
   plt.title('Phi = {:.1f} $\pi$'.format(num / np.pi))
   plt.ylabel('z/α')
   #plt.colorbar()
   return scat,

anim2 = animation.FuncAnimation(fig, animate2, frames = num_range, fargs=(scat, r, theta),
                              interval=500, blit=False)

# Circle1 = plt.Circle((0, 0), 1, color = 'red', fill = False)
# ax.add_patch(Circle1)

plt.show()
#%%
anim2.save('Harry/colormap_test.gif')
#%%

xn, yn, zn, un, vn, wn = loadBField('Output/dipole_nick_fieldaligned.npz')

xc, yc, zc, uc, vc, wc = reshapeCN(x, y, z, u, v, w, 4, 4, 4)

#%%

print(un - uc)

#print(x, y, z, u, v, w)