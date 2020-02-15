# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 15:06:22 2020

@author: Alexandre
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def FDTD_eq(E,H,E_temp,H_temp,t_index,Rb,Ra,Ca_arr,Cb_arr):
    # 3.93 a)
    # bad values are [:,-1,-1]
    H_temp[:,:,:,0,t_index] = H[:,:,:,0] + Rb*(np.roll(E[:,:,:,1],-1,axis=2) - E[:,:,:,1] - np.roll(E[:,:,:,2],-1,axis=1) + E[:,:,:,2])
    #H_temp[:,:,:,0] = H[:,:,:,0] + Rb*(np.roll(E[:,:,:,2],-1,axis=1) + E[:,:,:,2])

    # 3.93 b)
    # bad values are [-1,:,-1]
    H_temp[:,:,:,1,t_index] = H[:,:,:,1] + Rb*(np.roll(E[:,:,:,2],-1,axis=0) - E[:,:,:,2] - np.roll(E[:,:,:,0],-1,axis=2) + E[:,:,:,0])
    #H_temp[:,:,:,1] = H[:,:,:,1] + Rb*(np.roll(E[:,:,:,2],-1,axis=0) - E[:,:,:,2])

    # 3.93 c)
    # bad values are [-1,-1,:]
    H_temp[:,:,:,2,t_index] = H[:,:,:,2] + Rb*(np.roll(E[:,:,:,0],-1,axis=1) - E[:,:,:,0] - np.roll(E[:,:,:,1],-1,axis=0) + E[:,:,:,1])

    # 3.93 d)
    # bad values are [:,0,0]
    E_temp[:,:,:,0,t_index] = Ca_arr*E[:,:,:,0] + Cb_arr/Rb*(H[:,:,:,2] - np.roll(H[:,:,:,2],1,axis=1) - H[:,:,:,1] + np.roll(H[:,:,:,1],1,axis=2))
    # 3.93 e)
    # bad values are [0,:,0]
    E_temp[:,:,:,1,t_index] = Ca_arr*E[:,:,:,1] + Cb_arr/Rb*(H[:,:,:,0] - np.roll(H[:,:,:,0],1,axis=2) - H[:,:,:,2] + np.roll(H[:,:,:,2],1,axis=0))

    # 3.93 f)
    # bad values are [0,0,:]
    E_temp[:,:,:,2,t_index] = Ca_arr*E[:,:,:,2] + Cb_arr/Rb*(H[:,:,:,1] - np.roll(H[:,:,:,1],1,axis=0) - H[:,:,:,0] + np.roll(H[:,:,:,0],1,axis=1))


    index1 = np.array([0,-1])
    index2 = np.array([1,-2])
    temp = H_temp[...,t_index]
    temp[index1,:,:,:] = H[index2,:,:,:]
    temp[:,index1,:,:] = H[:,index2,:,:]
    temp[:,:,index1,:] = H[:,:,index2,:]
    H_temp[...,t_index] = temp

    temp = E_temp[...,t_index]
    temp[index1,:,:,:] = E[index2,:,:,:]
    temp[:,index1,:,:] = E[:,index2,:,:]
    temp[:,:,index1,:] = E[:,:,index2,:]
    E_temp[...,t_index] = temp

    return E_temp, H_temp

Npix = 30
Nmax = 20000

#conductivity
sigma= [0,0]
#permitivity ratio
ep_r = [1,1]
c = 3*10**8
ep_o = 10**(-9)/(36*np.pi)
mu_o = 4*np.pi*10**(-7)
freq = 2.5*10**9
umax = c/np.sqrt(ep_r[1])
wavelength = umax/freq
dx = wavelength/20
dy = wavelength/20
dz = wavelength/20
dt = dx/(300*c)

R = dt/(2*ep_o)
Ra = (dt/dx)**2/(ep_o*mu_o)

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
radius = dx*Npix//2
x = dx*np.array(range(Npix))
y = dy*np.array(range(Npix))
z = dz*np.array(range(Npix))

media = np.zeros((Npix,Npix,Npix))

arr=np.array(np.meshgrid(x,y,z))
mask = ((arr[0]-cx)**2 + (arr[1]-cy)**2 + (arr[2]-cz)**2 ) < radius**2

media[mask] = 1

Ca_val = []
Cb_val = []
Ca_val.append((1 - (R*sigma[0]/ep_r[0]))/(1 + (R*sigma[0]/ep_r[0])))
Ca_val.append((1 - (R*sigma[1]/ep_r[1]))/(1 + (R*sigma[1]/ep_r[1])))
Cb_val.append(Ra/(ep_r[0] + R*sigma[0]))
Cb_val.append(Ra/(ep_r[1] + R*sigma[1]))

Cb_arr = np.zeros((Npix,Npix,Npix))
Ca_arr = np.zeros((Npix,Npix,Npix))

Ca_arr[media == 0] = Ca_val[0]
Ca_arr[media == 1] = Ca_val[1]
Cb_arr[media == 0] = Cb_val[0]
Cb_arr[media == 1] = Cb_val[1]

E = np.zeros((Npix,Npix,Npix,3))
H = np.zeros((Npix,Npix,Npix,3))
print("Ca values : ",Ca_val)
print("Cb values : ",Cb_val)

num_save_steps = 3

E_temp = np.zeros((Npix,Npix,Npix,3,num_save_steps))
H_temp = np.zeros((Npix,Npix,Npix,3,num_save_steps))

# Store initial wave 
E[:,4,:,2] += np.cos(2*np.pi*freq*(0*dt))
E_temp[:,:,:,:,0] = E
H_temp[:,:,:,:,0] = H

t_index = 1
E_temp, H_temp = FDTD_eq(E,H,E_temp,H_temp,t_index,Rb,Ra,Ca_arr,Cb_arr)

t_index = 2

slices = np.zeros((Npix,Npix,Nmax//100))

for n in range(Nmax):
    print(n)
 
    E_temp, H_temp = FDTD_eq(E,H,E_temp,H_temp,t_index,Rb,Ra,Ca_arr,Cb_arr)

    # Set the n time step equal to the n-1 time step
    E = E_temp[:,:,:,:,1]
    H = H_temp[:,:,:,:,1]

    E_temp[:,:,:,:,1] = E_temp[:,:,:,:,2]
    H_temp[:,:,:,:,1] = H_temp[:,:,:,:,2]

    E[:,4,:,2] = np.cos(2*np.pi*freq*(n*dt))

    if n % 100 == 0:
        z = E[:,:,Npix//5,2]
        slices[:,:,n//100] = z

import matplotlib.animation as animation
from mpl_toolkits import mplot3d
"""
fig, ax = plt.subplots(projection='3d')

line, = ax.plot_surface(slices[:,:,0])

def animate(i):
    line.set_ydata(slices[:,:,i])
    #ax.set_ylim(-3,3)  # update the data.
    return line,


ani = animation.FuncAnimation(
    fig, animate, interval=40, blit=False, save_count=1000)

ani.save("movie.mp4")
"""

def update_plot(frame_number, zarray, plot):
    plot[0].remove()
    plot[0] = ax.plot_surface(x, y, slices[:,:,frame_number],cmap='seismic',linewidth=10)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x, y = np.meshgrid(x, y)

plot = [ax.plot_surface(x, y, slices[:,:,0], color='0.75', rstride=1, cstride=1,linewidth=10)]
ax.set_zlim(-3,3)
ani = animation.FuncAnimation(fig, update_plot,fargs=(slices, plot), interval=60,save_count=1000)
ani.save("movie.mp4")