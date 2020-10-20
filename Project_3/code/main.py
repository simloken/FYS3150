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
        self.a = 0
        self.a2 = 0
        
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

    def Verlet(self, primary, N, dt, beta=2, rel=False):
        v = np.zeros((N,2)); p = np.zeros((N,2))
        v[0] = self.v0; p[0] = self.r0
        dt2 = 0.5*dt*dt #save FLOPs
        if rel == False:
            Fx,Fy = self.force(primary, beta) #calculating the first acceleration here
        else:
            Fx,Fy = self.ForceRelativistic(primary,beta, v[0]) #relativistic force for i
        ax, ay = Fx/self.mass, Fy/self.mass #this saves us upwards of 50%+ runtime
        a = np.array((ax,ay)) #as we don't have to calculate acceleration twice in-loop
        for i in range(N-1):
            p[i+1] = p[i] + v[i]*dt+a*dt2
            self.r0 = p[i+1]
            if rel == False:
                Fx,Fy = self.force(primary,beta)
            else:
                Fx,Fy = self.ForceRelativistic(primary,beta, v[i])
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

    def ForceMultiBody(self, bodies, beta): #finds the x-and y-force on self exerted by primary
        TFx = 0; TFy = 0
        for primary in bodies:
            if self == primary:
                continue
            sx, sy = self.r0[0], self.r0[1] #self position  
            px, py = primary.r0[0], primary.r0[1] #primary position
            rx = px-sx
            ry = py-sy
            r = np.sqrt(rx**2+ry**2) #find radius from planet and it's primary
            F = G*self.mass*primary.mass/(r**beta) #force between planet and it's primary
            theta = atan2(ry,rx)
            Fx = F*np.cos(theta); Fy = F*np.sin(theta)
            TFx += Fx; TFy += Fy
        return TFx, TFy
    
    def VerletMultiBody(bodies, N, dt, beta=2): #not DONE
        dt2 = 0.5*dt*dt #save FLOPs
        for body in bodies:
            body.v = np.zeros((N,2)); body.p = np.zeros((N,2))
            body.v[0] = body.v0; body.p[0] = body.r0
            Fx,Fy = body.ForceMultiBody(bodies, beta)
            ax, ay = Fx/body.mass, Fy/body.mass
            body.a = np.array((ax,ay))
        for i in range(N-1):
            for body in bodies:
                body.p[i+1] = body.p[i] + body.v[i]*dt+body.a*dt2
                body.r0 = body.p[i+1]
                Fx,Fy = body.ForceMultiBody(bodies,beta)
                ax, ay = Fx/body.mass, Fy/body.mass
                body.a2 = body.a #old a
                body.a = np.array((ax,ay)) #forward a
                body.v[i+1] = body.v[i] + dt*((body.a+body.a2)/2)
            body.method = 'Verlet'

    def ForceRelativistic(self, primary, beta, vVec): #finds the x-and y-force on self exerted by primary
        
        sx, sy = self.r0[0], self.r0[1] #self position  
        px, py = primary.r0[0], primary.r0[1] #primary position
        rx = px-sx
        ry = py-sy
        rv = np.array((rx,ry))
        crs = np.cross(rv,vVec)
        c = 3e8
        r = np.sqrt(rx**2+ry**2) #find radius from planet and it's primary
        F = (G*self.mass*primary.mass/(r**beta))*(1+((3*crs**2)/(r**2*c**2))) #force between planet and it's primary
        theta = atan2(ry,rx)
        Fx = F*np.cos(theta); Fy = F*np.sin(theta)
        return Fx, Fy
    
    def VerletPerihelion(self, primary, dt, beta=2, rel=True):
        N = 1000
        v = np.zeros((N,2)); p = np.zeros((N,2))
        v[0] = self.v0; p[0] = self.r0
        dt2 = 0.5*dt*dt #save FLOPs
        if rel == False:
            Fx,Fy = self.force(primary, beta) #calculating the first acceleration here
        else:
            Fx,Fy = self.ForceRelativistic(primary,beta, v[0]) #relativistic force for i
        ax, ay = Fx/self.mass, Fy/self.mass #this saves us upwards of 50%+ runtime
        a = np.array((ax,ay)) #as we don't have to calculate acceleration twice in-loop
        i = 0
        j = 0
        k = 0
        while k == 0:
            p[i+1] = p[i] + v[i]*dt+a*dt2
            p1 = np.linalg.norm(p[i])
            p2 = np.linalg.norm(p[i+1])
            
            if p2 < p1 and j == 0: #sinking
                j = 1
                print('Sinking')
            if j == 1:
                if p2 > p1: #increasing
                    print('Solution found')
                    return p[i]
                    k = 1
                    
            self.r0 = p[i+1]
            if rel == False:
                Fx,Fy = self.force(primary,beta)
            else:
                Fx,Fy = self.ForceRelativistic(primary,beta, v[i])
            ax, ay = Fx/self.mass, Fy/self.mass
            a2 = a #old a
            a = np.array((ax,ay)) #forward a
            v[i+1] = v[i] + dt*((a+a2)/2)
            i += 1
        self.p = p
        self.v = v
        self.method = 'Verlet'
#b - DONE    
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

#e - NOT DONE?
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



#g - DONE
earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,30e3)), 6e24) #re-initializing
sun = CelestialBody('Sun', np.array((0,0)), np.array((0,0)), 2e30) #re-initializing
jupiter = CelestialBody('Jupiter', np.array((0,5.2*AU)), np.array((13e3,0)), 1.9e27)
bodies= [earth, jupiter, sun]
CelestialBody.VerletMultiBody(bodies, 30*365, 24*3600)
plt.plot(earth.p[:,0]/AU,earth.p[:,1]/AU, label='Earth')
plt.plot(jupiter.p[:,0]/AU, jupiter.p[:,1]/AU, label='Jupiter')
plt.plot(sun.p[:,0]/AU, sun.p[:,1]/AU, label='Sun')
plt.xlabel('x position [AU]')
plt.ylabel('y position [AU]')
plt.title('Earth and Jupiters Orbit with real Jupiter mass over 30 years')
plt.legend()
plt.show()

#h
comx = (1.9e27*5.2*AU)/(2e30+6e24+1.9e27)
comy = (6e24*AU)/(2e30+6e24+1.9e27)
sun = CelestialBody('Sun', np.array((-comx,-comy)), np.array((-12.35,-0.09)), 2e30) #re-initializing
earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,30e3)), 6e24) #re-initializing
jupiter = CelestialBody('Jupiter', np.array((0,5.2*AU)), np.array((13e3,0)), 1.9e27)
bodies = [earth,jupiter, sun]
CelestialBody.VerletMultiBody(bodies, 50*365, 24*3600)
plt.plot(earth.p[:,0]/AU,earth.p[:,1]/AU, label='Earth')
plt.plot(jupiter.p[:,0]/AU, jupiter.p[:,1]/AU, label='Jupiter')
plt.plot(sun.p[:,0]/AU, sun.p[:,1]/AU, label='Sun')
plt.xlabel('x position [AU]')
plt.ylabel('y position [AU]')
plt.title('Earth and Jupiters Orbit with real Jupiter mass over 50 years\nwith initialized sun with regards to Center of Mass')
plt.legend()
plt.show()

""" center of mass for system with all 8 planets (and Pluto) (antiquated, but kept around as writing out all of this takes a while)
comx = (AU*(-3.51e-1*3.3e23+-2.4e-1*4.9e24+9.05e-1*6e24+1.3*6.6e23+2.56*1.9e27+5.15*5.5e26+1.55e1*8.8e25+2.94e1*1.03e26+1.38e1*1.3e22)/
        (2e30+3.3e23+4.9e24+6e24+6.6e23+1.9e27+5.5e26+8.8e25+1.03e26+1.3e22))
comy = (AU*(-6.78e-2*3.3e23+6.85e-1*4.9e24+4.09e-1*6e24+5.52e-1*6.6e23-4.42*1.9e27-8.56*5.5e26+1.22e1*8.8e25-5.46*1.03e26-3.13e1*1.3e22)/
        (2e30+3.3e23+4.9e24+6e24+6.6e23+1.9e27+5.5e26+8.8e25+1.03e26+1.3e22))"""

sun = CelestialBody('Sun', AU*np.array((-6.11e-3,6.42e-3)), np.array((-0.02e3,-0.009e3)), 2e30) #re-initializing
mercury = CelestialBody('Mercury', AU*np.array((3.54e-1,-6.78e-2)), np.array((0.38e3,50e3)), 3.3e23) #re-initializing
venus = CelestialBody('Venus', AU*np.array((-2.4e-1,6.85e-1)), np.array((-33e3,-11.6e3)), 4.9e24) #re-initializing
earth = CelestialBody('Earth', AU*np.array((9.05e-1,4.09e-1)), np.array((-12.5e3,27.2e3)), 6e24) #re-initializing
mars = CelestialBody('Mars', AU*np.array((1.3,5.52e-1)), np.array((-8.4e3,24.4e3)), 6.6e23) #re-initializing
jupiter = CelestialBody('Jupiter', AU*np.array((2.56,-4.42)), np.array((11.1e3,7.17e3)), 1.9e27)
saturn = CelestialBody('Saturn', AU*np.array((5.15,-8.56)), np.array((7.74e3,4.95e3)), 5.5e26) #re-initializing
uranus = CelestialBody('Uranus', AU*np.array((1.55e1,1.22e1)), np.array((-4.26e3,5e3)), 8.8e25) #re-initializing
neptune = CelestialBody('Neptune', AU*np.array((2.94e1,-5.46)), np.array((0.96e3,5.38e3)), 1.03e26) #re-initializing
#pluto = CelestialBody('Pluto', AU*np.array((1.38e1,-3.12e1)), np.array((5.11e3,1.04e3)), 1.3e22) #re-initializing
bodies = [sun,mercury,venus,earth,mars,jupiter,saturn,uranus,neptune]
CelestialBody.VerletMultiBody(bodies, 200*365, 24*3600)
plt.figure(figsize=(10,10))
for body in bodies:
    plt.plot(body.p[:,0]/AU, body.p[:,1]/AU, label='%s' %(body.name))
plt.xlabel('x position [AU]')
plt.ylabel('y position [AU]')
plt.title('Solar System Orbits over 200 years')
plt.legend(loc='upper left')
plt.show()

#i - NOT DONE

sun = CelestialBody('Sun', np.array((0,0)), np.array((0,0)), 2e30) #re-initializing
mercury = CelestialBody('Mercury', AU*np.array((0,0.3075)), np.array((58.97e3,0)), 3.3e23) #re-initializing
mercury.Verlet(sun,100*365,24*3600, rel=True)
ps = mercury.p
plows = np.zeros(len(ps))
for i in range(len(ps)):
    plows[i] = np.linalg.norm(ps[i])
mini = np.argmin(plows[1:100])
thetap1 = np.arctan(ps[mini,1]/ps[mini,0])
mercury = CelestialBody('Mercury', np.array((mercury.p[100*365-1,0],mercury.p[100*365-1,1])),
                        np.array((mercury.v[100*365-1,0],mercury.v[100*365-1,1])), 3.3e23) #re-initializing
century_p = mercury.VerletPerihelion(sun,24*3600)
thetap2 = np.arctan(century_p[1]/century_p[0])
print(thetap1, thetap2)
