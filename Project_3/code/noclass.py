# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from math import atan2
M_earth = 6e24
M_sun = 2e30
r0 = 149.6e9  #inital pos
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

plt.plot(x[:,0],x[:,1], label='Earth')
plt.legend()
plt.show()