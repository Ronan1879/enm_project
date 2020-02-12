# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 15:06:22 2020

@author: Alexandre
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

Npix = 100
Nmax = 100

#conductivity
sigma= [0,0]
#permitivity ratio
ep_r = [1,2]
c = 3*10**8
ep_o = 8.85*10**(-12)
mu_o = 12.57*10**(-7)
freq = 2.5*10**9
umax = c/np.sqrt(ep_r[1])
wavelength = umax/freq
dx = dy = dz = wavelength/10
dt = dx/(c)

R = dt/(2*ep_o)
Ra = (c*dt/(dx))**2

Rb = dt/(mu_o*dx)

print("wavelength : ",wavelength)
print("umax : ",umax)
print("spatial increment : ",dx)
print("time increment : ",dt)
print("R : ",R)
print("Ra : ",Ra)
print("Rb : ",Rb)

print("Increment condition respected? : ",umax*dt/dx < 1/np.sqrt(3),umax*dt/dx)

# Define x and y array and find center of grid cx, cy
cx = dx*Npix//2
cy = dy*Npix//2
cz = dz*Npix//2
radius = dx*Npix//4
x = dx*np.array(range(Npix))
y = dy*np.array(range(Npix))
z = dz*np.array(range(Npix))

media = np.zeros((Npix,Npix,Npix))

arr=np.array(np.meshgrid(x,y,z))
mask = ((arr[0]-cx)**2 + (arr[1]-cy)**2 + (arr[2]-cz)**2 ) <= radius**2


media[mask] = 1

slices = []
Ca_val = []
Cb_val = []
Ca_val.append((1 - (R*sigma[0]/ep_r[0]))/(1 + (R*sigma[0]/ep_r[0])))
Ca_val.append((1 - (R*sigma[1]/ep_r[1]))/(1 + (R*sigma[1]/ep_r[1])))
Cb_val.append(Ra/(ep_r[0] + R*sigma[0]))
Cb_val.append(Ra/(ep_r[1] + R*sigma[1]))

Cb_arr = Ca_arr = np.zeros((Npix,Npix,Npix))

Ca_arr[media == 0] = Ca_val[0]
Ca_arr[media == 1] = Ca_val[1]
Cb_arr[media == 0] = Cb_val[0]
Cb_arr[media == 1] = Cb_val[1]

E = H = np.zeros((Npix,Npix,Npix,3))

print("Ca values : ",Ca_val)
print("Cb values : ",Cb_val)

for n in range(Nmax):
    print(n)
    E_temp = H_temp = np.zeros((Npix,Npix,Npix,3))
    
    E[:,5,:,2] += np.cos(2*np.pi*freq*(n*dt))

    # 3.93 a)
    # bad values are [:,-1,-1]
    H_temp[:,:,:,0] = H[:,:,:,0] + Rb*(np.roll(E[:,:,:,1],-1,axis=2) - E[:,:,:,1] - np.roll(E[:,:,:,2],-1,axis=1) + E[:,:,:,2])

    # 3.93 b)
    # bad values are [-1,:,-1]
    H_temp[:,:,:,1] = H[:,:,:,1] + Rb*(np.roll(E[:,:,:,2],-1,axis=0) - E[:,:,:,2] - np.roll(E[:,:,:,0],-1,axis=2) + E[:,:,:,0])

    # 3.93 c)
    # bad values are [-1,-1,:]
    H_temp[:,:,:,2] = H[:,:,:,2] + Rb*(np.roll(E[:,:,:,0],-1,axis=1) - E[:,:,:,0] - np.roll(E[:,:,:,1],-1,axis=0) + E[:,:,:,1])

    # 3.93 d)
    # bad values are [:,0,0]
    E_temp[:,:,:,0] = Ca_arr*Rb*E[:,:,:,0] + Cb_arr*(H[:,:,:,2] - np.roll(H[:,:,:,2],1,axis=1) - H[:,:,:,1] + np.roll(H[:,:,:,1],1,axis=2))
    # 3.93 e)
    # bad values are [0,:,0]
    E_temp[:,:,:,1] = Ca_arr*Rb*E[:,:,:,1] + Cb_arr*(H[:,:,:,0] - np.roll(H[:,:,:,0],1,axis=2) - H[:,:,:,2] + np.roll(H[:,:,:,2],1,axis=0))

    # 3.93 f)
    # bad values are [0,0,:]
    E_temp[:,:,:,2] = Ca_arr*Rb*E[:,:,:,2] + Cb_arr*(H[:,:,:,1] - np.roll(H[:,:,:,1],1,axis=0) - H[:,:,:,0] + np.roll(H[:,:,:,0],1,axis=1))

    E_temp *= Rb
    """
    index1 = np.array([0,-1])
    index2 = np.array([1,-2])
    H_temp[index1,:,:] = H[index2,:,:]
    E_temp[index1,:,:] = E[index2,:,:]
    H_temp[:,index1,:] = H[:,index2,:]
    E_temp[:,index1,:] = E[:,index2,:]
    H_temp[:,:,index1] = H[:,:,index2]
    E_temp[:,:,index1] = E[:,:,index2]
    """
    E=E_temp
    H=H_temp
    

    z = E[Npix//4,:,Npix//5,2]
    slices.append(z)

import matplotlib.animation as animation

fig, ax = plt.subplots()

line, = ax.plot(slices[0])

def animate(i):
    line.set_ydata(slices[i])
    ax.set_ylim(-2*np.max(slices[0]),2*np.max(slices[0]))  # update the data.
    return line,


ani = animation.FuncAnimation(
    fig, animate, interval=50, blit=False, save_count=500)

ani.save("movie.mp4")
