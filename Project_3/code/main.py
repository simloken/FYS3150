# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from math import atan2
class CelestialBody:
    def __init__(self, name, r0, v0, mass):
        
        if isinstance(name, str) == False: #ensure name is a string
            raise TypeError('name must be str')
        
        if hasattr(r0, "__len__") == True: #strings have len, ensure r0 is not str
            if isinstance(r0, str) == True:
                raise TypeError('r0 must be array-like, not string')
            elif r0.any() == str:
                raise TypeError('r0 must be array-like, not string')
        if hasattr(v0, "__len__") == True: #strings have len, ensure v0 is not str
            if isinstance(v0, str) == True:
                raise TypeError('v0 must be array-like, not string')
            elif v0.any() == str:
                raise TypeError('v0 must be array-like, not string')
        if hasattr(r0, "__len__") == False:  #ensure r0 is list or array
            raise ValueError('r0 must be array-like')
        if hasattr(v0, "__len__") == False: #ensure v0 is list or array
            raise ValueError('v0 must be array-like')
        if len(r0) > 2: #ensure r0 is shape (x,y)
            raise ValueError('r0 can have no more than two elements')
        if len(np.shape(r0)) != 1: #ensure r0 is shape (x,y)
            raise ValueError('r0 can have no more than two elements')
        if len(v0) > 2: #ensure v0 is shape (x,y)
            raise ValueError('v0 can have no more than two elements')   
        if len(np.shape(v0)) != 1: #ensure v0 is shape (x,y)
            raise ValueError('v0 can have no more than two elements')
        
        if isinstance(mass, int) == False and isinstance(mass, float) == False: #ensure mass is an integer or a float
            raise TypeError('mass must be float or int')
            
        self.name = name
        self.r0 = r0
        self.v0 = v0
        self.mass = mass
        self.p = 0
        self.v = 0
        
    def force(self, primary): #finds the x-and y-force on self exerted by primary
        
        sx, sy = self.r0[0], self.r0[1] #self position  
        px, py = primary.r0[0], primary.r0[1] #primary position
        rx = px-sx
        ry = py-sy
        r = np.sqrt(rx**2+ry**2) #find radius from planet and it's primary
        F = G*self.mass*primary.mass/(r**2) #force between planet and it's primary
        theta = atan2(ry,rx)
        Fx = F*np.cos(theta); Fy = F*np.sin(theta)
        return Fx, Fy
        
    
    def fEuler(self, primary, N,dt):
        v = np.zeros((N,2)); p = np.zeros((N,2))
        v[0] = self.v0; p[0] = self.r0
        for i in range(N-1):
            Fx,Fy = self.force(primary)
            ax, ay = Fx/self.mass, Fy/self.mass
            a = np.array((ax,ay))
            v[i+1] = v[i] + a*dt
            p[i+1] = p[i] + v[i+1]*dt
            self.r0 = p[i+1]
        self.p = p
        self.v = v

    def Verlet():
        lol = 'lmao'
        
G = 6.67e-11
AU = 149.6e9       
earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,30e3)), 6e24)
sun = CelestialBody('Sun', np.array((0,0)), np.array((0,0)), 2e30)

earth.fEuler(sun,365,24*3600)
l = earth.p
plt.plot(earth.p[:,0],earth.p[:,1], label='%s' %(earth.name))
plt.legend()
plt.show()