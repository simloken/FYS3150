# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from math import atan2
def proj_3a():
    M_earth = 6e24
    M_sun = 2e30
    AU = 149.6e9
    r0 = AU  #inital pos
    v0 = 30e3
    G = 6.67e-11
    N = 365
    dt = 24*3600
    x = np.zeros((N,2))
    v = np.zeros((N,2))
    x[0,0] = r0
    v[0,1] = v0
    for i in range(N-1):
        r = np.sqrt(x[i,0]**2+x[i,1]**2)
        FG = G*M_sun*M_earth/(r**2)
        theta = atan2(x[i,1],x[i,0])
        fx = -FG*np.cos(theta)
        fy = -FG*np.sin(theta)
        ax = fx/M_earth; ay = fy/M_earth
        a = np.array((ax,ay))
        v[i+1] = v[i] + a*dt
        x[i+1] = x[i] + v[i+1]*dt
    
    
    v2 = np.zeros((N,2)); x2 = np.zeros((N,2))
    v2[0,1] = v0; x2[0,0] = r0
    dt2 = 0.5*dt*dt #save FLOPs
    r = np.linalg.norm(x2[0])
    FG = G*M_sun*M_earth/(r**2)
    theta = atan2(x2[0,1],x2[0,0])
    fx = -FG*np.cos(theta)
    fy = -FG*np.sin(theta)
    ax, ay = fx/M_earth, fy/M_earth
    a = np.array((ax,ay))
    for i in range(N-1):
        x2[i+1] = x2[i] + v2[i]*dt+a*dt2
        r = np.linalg.norm(x2[i+1])
        FG = G*M_sun*M_earth/(r**2)
        theta = atan2(x2[i+1,1],x2[i+1,0])
        fx = -FG*np.cos(theta)
        fy = -FG*np.sin(theta)
        ax = fx/M_earth; ay = fy/M_earth
        a2 = a #old a
        a = np.array((ax,ay)) #forward a
        v2[i+1] = v2[i] + dt*((a+a2)/2)
    
    plt.plot(x[:,0]/AU,x[:,1]/AU, label='Earth w/ Forward Euler')
    plt.plot(x2[:,0]/AU,x2[:,1]/AU, label='Earth w/ Verlet Method')
    plt.legend()
    plt.xlabel('x position [AU]')
    plt.ylabel('y position [AU]')
    plt.show()