# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from math import atan2
import time
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
        self.method = 0
        
    def force(self, primary, beta): #finds the x-and y-force on self exerted by primary
        
        sx, sy = self.r0[0], self.r0[1] #self position  
        px, py = primary.r0[0], primary.r0[1] #primary position
        rx = px-sx
        ry = py-sy
        r = np.sqrt(rx**2+ry**2) #find radius from planet and it's primary
        F = G*self.mass*primary.mass/(r**beta) #force between planet and it's primary
        theta = atan2(ry,rx)
        Fx = F*np.cos(theta); Fy = F*np.sin(theta)
        return Fx, Fy
        
    
    def fEuler(self, primary, N,dt):
        beta = 2
        v = np.zeros((N,2)); p = np.zeros((N,2))
        v[0] = self.v0; p[0] = self.r0
        for i in range(N-1):
            Fx,Fy = self.force(primary, beta)
            ax, ay = Fx/self.mass, Fy/self.mass
            a = np.array((ax,ay))
            v[i+1] = v[i] + a*dt
            p[i+1] = p[i] + v[i+1]*dt
            self.r0 = p[i+1]
        self.p = p
        self.v = v
        self.method = 'Forward Euler'

    def Verlet(self, primary, N, dt, beta=2):
        v = np.zeros((N,2)); p = np.zeros((N,2))
        v[0] = self.v0; p[0] = self.r0
        dt2 = 0.5*dt*dt #save FLOPs
        Fx,Fy = self.force(primary, beta) #calculating the first acceleration here
        ax, ay = Fx/self.mass, Fy/self.mass #this saves us upwards of 50%+ runtime
        a = np.array((ax,ay)) #as we don't have to calculate acceleration twice in-loop
        for i in range(N-1):
            p[i+1] = p[i] + v[i]*dt+a*dt2
            self.r0 = p[i+1]
            Fx,Fy = self.force(primary,beta)
            ax, ay = Fx/self.mass, Fy/self.mass
            a2 = a #old a
            a = np.array((ax,ay)) #forward a
            v[i+1] = v[i] + dt*((a+a2)/2)
        self.p = p
        self.v = v
        self.method = 'Verlet'
      
    def Energies(self,primary,N,dt):
        sx, sy = self.r0[0], self.r0[1] #self position  
        px, py = primary.r0[0], primary.r0[1] #primary position
        rx = px-sx
        ry = py-sy
        r = np.sqrt(rx**2+ry**2) #find radius from planet and it's primary
        
        ETOT = np.zeros(N); U = np.zeros(N); K = np.zeros(N)
        U[0] = -G*self.mass*primary.mass/r
        K[0] = 0.5*self.mass*np.linalg.norm(self.v0)**2
        ETOT[0] = U[0]+K[0]
        
        v = np.zeros((N,2)); p = np.zeros((N,2))
        v[0] = self.v0; p[0] = self.r0
        beta = 2
        dt2 = 0.5*dt*dt #save FLOPs
        Fx,Fy = self.force(primary, beta) #calculating the first acceleration here
        ax, ay = Fx/self.mass, Fy/self.mass #this saves us upwards of 50%+ runtime
        a = np.array((ax,ay)) #as we don't have to calculate acceleration twice in-loop
        for i in range(N-1):
            p[i+1] = p[i] + v[i]*dt+a*dt2
            self.r0 = p[i+1]
            Fx,Fy = self.force(primary, beta)
            ax, ay = Fx/self.mass, Fy/self.mass
            a2 = a #old a
            a = np.array((ax,ay)) #forward a
            v[i+1] = v[i] + dt*((a+a2)/2)
            U[i+1] = -G*self.mass*primary.mass/np.linalg.norm(self.r0)
            K[i+1] = 0.5*self.mass*np.linalg.norm(v[i+1])**2
            ETOT[i+1] = U[i+1]+K[i+1]
        return U, K, ETOT
"""   
    def VerletMultiBody(self, primary, N, dt, beta=2): #not DONE
        v = np.zeros((N,2)); p = np.zeros((N,2))
        v[0] = self.v0; p[0] = self.r0
        dt2 = 0.5*dt*dt #save FLOPs
        TFx, TFy = 0
        for bodies in something:
            Fx,Fy = self.force(bodies, beta) #calculating the first acceleration 
            TFx, TFy += Fx, Fy
        ax, ay = Fx/self.mass, Fy/self.mass #this saves us upwards of 50%+ runtime
        a = np.array((ax,ay)) #as we don't have to calculate acceleration twice in-loop
        for i in range(N-1):
            p[i+1] = p[i] + v[i]*dt+a*dt2
            self.r0 = p[i+1]
            TFx, TFy = 0
            for bodies in something:
                Fx,Fy = self.force(bodies,beta)
                TFx, TFy += Fx, Fy
            ax, ay = Fx/self.mass, Fy/self.mass
            a2 = a #old a
            a = np.array((ax,ay)) #forward a
            v[i+1] = v[i] + dt*((a+a2)/2)
        self.p = p
        self.v = v
        self.method = 'Verlet'
"""  
#B - DONE    
G = 6.67e-11; AU = 149.6e9
earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,30e3)), 6e24)
sun = CelestialBody('Sun', np.array((0,0)), np.array((0,0)), 2e30)

t = time.perf_counter()
earth.fEuler(sun,365,24*3600)
t1 = time.perf_counter() -t
E = earth; E.p = E.p/AU #scaling

earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,30e3)), 6e24) #re-initializing

t = time.perf_counter()
earth.Verlet(sun,365,24*3600)
t2 = time.perf_counter() - t
V = earth; V.p = V.p/AU #scaling


print('Time elapsed for Forward Euler Method:', t1)
print('Time elapsed for Verlet Method:', t2)



plt.plot(E.p[:,0],E.p[:,1], label='%s using %s method' %(E.name, E.method))
plt.plot(V.p[:,0],V.p[:,1], label='%s using %s method' %(V.name, V.method))
plt.xlabel('x position [AU]')
plt.ylabel('y position [AU]')
plt.title('Earths orbit around the sun')
plt.legend()
plt.show()
#c - DONE
earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,30e3)), 6e24) #re-initializing
U, K, ETOT = earth.Energies(sun,365,24*3600)

plt.figure()
plt.title('All energies of %s in a complete orbit' %(earth.name))
plt.plot(range(len(ETOT)),U, label='Potential Energy')
plt.plot(range(len(ETOT)),K, label='Kinetic Energy')
plt.plot(range(len(ETOT)),ETOT, label='Total Energy')
plt.xlabel('Days')
plt.ylabel('Energy [J]')
plt.legend()
plt.show()
plt.figure()
plt.title('Total Energy of %s in a complete orbit' %(earth.name))
plt.plot(range(len(ETOT)),ETOT, label='Total Energy')
plt.xlabel('Days')
plt.ylabel('Energy [J]')
plt.legend()
plt.show()

#e - NOT DONE
betas = np.linspace(2,3,11)
years = 10
bOrbits = np.zeros((len(betas),years*365,2))
j = 0
for i in betas:
    earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,30e3)), 6e24) #re-initializing for each loop
    earth.Verlet(sun,years*365,48*3600,i)
    bOrbits[j] = earth.p
    j += 1
bOrbits = bOrbits/AU
plt.figure()
for i in range(len(betas)):
    b = betas[i]
    plt.plot(bOrbits[i,:,0], bOrbits[i,:,1], label='B = %.1f' %(b))
plt.legend()
plt.show()


#f - DONE
expected_escape_vel = np.sqrt(2*G*sun.mass/AU)
print(expected_escape_vel)
vels = np.linspace(35e3,45e3,11)
plt.figure(figsize=(7.2,7.2))
for i in vels:
    earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,i)), 6e24) #re-initializing for each loop
    earth.Verlet(sun,50*365,24*3600) #orbit over 50 years
    plt.plot(earth.p[:,0]/AU,earth.p[:,1]/AU, label='Orbit for intial velocity %.0f m/s' %(i))
plt.legend()
plt.xlabel('x position [AU]')
plt.ylabel('y position [AU]')
plt.title('Earths orbit around the sun for different initial velocities')
plt.show()



#g - NOT DONE
earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,30e3)), 6e24) #re-initializing
sun = CelestialBody('Sun', np.array((0,0)), np.array((0,0)), 2e30) #re-initializing
jupiter = CelestialBody('Jupiter', np.array((0,5.2*AU)), np.array((13e3,0)), 1.9e27)
