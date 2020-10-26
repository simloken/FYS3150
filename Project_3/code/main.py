# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from noclass import proj_3a
import time
"""
The main class used in essentially all of project 3.
Initializes given initial coordinates, initial velocity, mass and a name for our object
Is only compatible with two dimensions, though this could easily be modified.
Ensures no wrong inputs are given, like too large an array for position, and
that all values are the correct type.
"""
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
    
    """
    Basic force calculation using the supplied graviational force equation.
    Takes three inputs, the two bodies the force is acting between and a beta
    value which is passed from it's corresponding solving method (Verlet/Euler).
    """    
    def force(self, primary, beta): #finds the x-and y-force on self exerted by primary
        G = 6.67e-11
        sx, sy = self.r0[0], self.r0[1] #self position  
        px, py = primary.r0[0], primary.r0[1] #primary position
        rx = px-sx
        ry = py-sy
        r = np.linalg.norm((rx,ry)) #find radius from planet and it's primary
        F = G*self.mass*primary.mass/(r**beta) #force between planet and it's primary
        theta = np.arctan2(ry,rx)
        Fx = F*np.cos(theta); Fy = F*np.sin(theta)
        return Fx, Fy
        
    """
    Takes the two bodies, amount of steps N and our timestep dt. This is also
    where beta is defined, which defaults to beta=2.
    Solves our system using the Forward Euler method
    When it's done, it assigns two arrays (position and velocity) to the body 
    in which we're plotting which we can then retrieve from our object to plot
    """
    def fEuler(self, primary, N,dt, beta=2):
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

    """
    Takes two bodies, amount of steps N and timestep dt.
    Additionally takes the previously discussed beta value and a boolean rel.
    Rel dictates whether or not we want to make relativistic corrections to our
    force function, and is only used in 3i
    Solves using the Verlet Velocity Method, and assigns two arrays onto our
    object which we can use to plot.
    """
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
            elif rel == True:
                Fx,Fy = self.ForceRelativistic(primary,beta, v[i])
            ax, ay = Fx/self.mass, Fy/self.mass
            a2 = a #old a
            a = np.array((ax,ay)) #forward a
            v[i+1] = v[i] + dt*((a+a2)/2)
        self.p = p
        self.v = v
        self.method = 'Verlet'
     
    """
    Takes two bodies, amount of steps N and timestep dt and the value beta
    Solves our system, same as before using the Verlet method.
    Does not assign arrays onto our object, instead returns three arrays:
    U: Potential Energy, K: Kinetic Energy and ETOT: Total Energy.
    """
    def Energies(self,primary,N,dt, beta=2):
        G = 6.67e-11
        sx, sy = self.r0[0], self.r0[1] #self position  
        px, py = primary.r0[0], primary.r0[1] #primary position
        rx = px-sx
        ry = py-sy
        r = np.linalg.norm((rx,ry)) #find radius from planet and it's primary
        ETOT = np.zeros(N); U = np.zeros(N); K = np.zeros(N)
        U[0] = -G*self.mass*primary.mass/r
        K[0] = 0.5*self.mass*np.linalg.norm(self.v0)**2
        ETOT[0] = U[0]+K[0]
        
        v = np.zeros((N,2)); p = np.zeros((N,2))
        v[0] = self.v0; p[0] = self.r0
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
    Takes two bodies, amount of steps N and timestep dt and the value beta
    Solves our system, same as before using the Verlet method.
    Does not assign arrays onto our object, instead returns an array l.
    l: Angular Momentum
    """
    def AngularMom(self,primary,N,dt,beta=2):
        v = np.zeros((N,2)); p = np.zeros((N,2)); l = np.zeros(N)
        v[0] = self.v0; p[0] = self.r0
        r = np.array((p[0,0], p[0,1])); mom = v[0]*self.mass
        l[0] = np.cross(r,mom)
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
            r = np.array((p[i+1,0], p[i+1,1])); mom = v[i+1]*self.mass
            l[i+1] = np.cross(r,mom)
        return l
    """
    Takes an array of all celestialbodies in our system and calculates gravitational
    force between them. Also takes a value beta.
    Makes sure not to calculate graviational force onto itself.
    """
    def ForceMultiBody(self, bodies, beta): #finds the x-and y-force on self exerted by primary
        G = 6.67e-11
        TFx = 0; TFy = 0
        for primary in bodies:
            if self == primary:
                continue
            sx, sy = self.r0[0], self.r0[1] #self position  
            px, py = primary.r0[0], primary.r0[1] #primary position
            rx = px-sx
            ry = py-sy
            r = np.linalg.norm((rx,ry)) #find radius from planet and it's primary
            F = G*self.mass*primary.mass/(r**beta) #force between planet and it's primary
            theta = np.arctan2(ry,rx)
            Fx = F*np.cos(theta); Fy = F*np.sin(theta)
            TFx += Fx; TFy += Fy
        return TFx, TFy
    """
    The corresponding solver for the multibody problem. Takes an array of all celestial bodies.
    Loops through them all and calculates all their movement using the Verlet Velocity method.
    The initial for-loop is used to create arrays, which are then filled in the
    second for-loop.
    Returns nothing, but assigns two arrays onto each celestial body which we
    can use to plot their orbits.
    """
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
    """
    Takes two acting bodies, a beta value and a velocity vector.
    This is the relativistic correction of the force which is called if and only
    if the boolean rel=True.
    Returns x-and y-forces respectively
    """
    def ForceRelativistic(self, primary, beta, vVec): #finds the x-and y-force on self exerted by primary
        G = 6.67e-11
        sx, sy = self.r0[0], self.r0[1] #self position  
        px, py = primary.r0[0], primary.r0[1] #primary position
        rx = px-sx
        ry = py-sy
        rv = np.array((rx,ry))
        l = abs(np.cross(rv,vVec))
        cc = 3e8**2
        r = np.linalg.norm((rx,ry)) #find radius from planet and it's primary
        F = (G*self.mass*primary.mass/(r**beta))*(1+((3*l**2)/(r**2*cc))) #force between planet and it's primary
        theta = np.arctan2(ry,rx)
        Fx = F*np.cos(theta); Fy = F*np.sin(theta)
        return Fx, Fy
    """
    A special verlet velocity method to be called AFTER one century has passed
    for Mercury's orbit.
    Keeps letting Mercury orbit until it's at the Perihelion, fromwhich it returns
    the coordinates.
    These can then be used to check if we have moved 42 arcseconds from our initial
    position when at the Perihelion
    """
    def VerletPerihelion(self, primary, dt, tol, beta=2, rel=True):
        AU = 149.6e9
        N = 100000
        v = np.zeros((N,2)); p = np.zeros((N,2))
        v[0] = self.v0; p[0] = self.r0
        dt2 = 0.5*dt*dt #save FLOPs
        if rel == False:
            Fx,Fy = self.force(primary, beta) #calculating the first acceleration here
        else:
            Fx,Fy = self.ForceRelativistic(primary,beta, v[0]) #relativistic force for i
        ax, ay = Fx/self.mass, Fy/self.mass #this saves us upwards of 50%+ runtime
        a = np.array((ax,ay)) #as we don't have to calculate acceleration twice in-loop
        pp = np.linalg.norm(p[0])
        if pp - 0.3075*AU < tol:
            self.p = p; self.v = v
            return p[0], 0
        oldp = pp
        i = 0
        j = 0
        k = 0
        while k == 0:
            p[i+1] = p[i] + v[i]*dt+a*dt2
            pp = np.linalg.norm(p[i+1])
            if j == 1 and pp > oldp: #next point is not closer than previous
                self.p = p; self.v = v
                return p[i+1], i
            else:
                j = 0
            if pp - 0.3075*AU < tol and j == 0: #within acceptable range
                        j = 1
            self.r0 = p[i+1]
            oldp = pp
            if rel == False:
                Fx,Fy = self.force(primary,beta)
            else:
                Fx,Fy = self.ForceRelativistic(primary,beta, v[i])
            ax, ay = Fx/self.mass, Fy/self.mass
            a2 = a #old a
            a = np.array((ax,ay)) #forward a
            v[i+1] = v[i] + dt*((a+a2)/2)
            i += 1
        
 
"""
Solves task 3b
Essentially the same as 3a, except now using OOP.
Additionally prints runtimes for the Forward Euler Method and Velocity Verlet Method respectively
Call the function from the console.
"""         
def proj_3b():
    AU = 149.6e9
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
    
    
    plt.figure()
    plt.plot(E.p[:,0],E.p[:,1], label='%s using %s method' %(E.name, E.method))
    plt.plot(V.p[:,0],V.p[:,1], label='%s using %s method' %(V.name, V.method))
    plt.xlabel('x position [AU]')
    plt.ylabel('y position [AU]')
    plt.title('Earths orbit around the Sun over one year\nusing two different methods (OOP)')
    plt.legend()
    plt.show()

"""
Solves task 3c
Checks the stability of our Velocity Verlet method for different timesteps
Plots the potential, kinetic and total energy of Earth over one year
Additionally plots a secondary plot with just the total energy by itself
Call the function from the console
"""
def proj_3c():
    AU = 149.6e9
    steps = [24,6,1,120]; j = 0
    fig, axs = plt.subplots(2,2)
    fig.set_figheight(10)
    fig.set_figwidth(10)
    sun = CelestialBody('Sun', np.array((0,0)), np.array((0,0)), 2e30)
    for i in steps:
        earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,30e3)), 6e24) #re-initializing
        earth.Verlet(sun,int(365*24/i),i*3600)
        p = earth.p
        if j == 0:
            axs[0,0].plot(p[:,0]/AU,p[:,1]/AU)
            axs[0,0].set_title('N = {}{} {}  $\Delta$t = {}{}' .format(int(365*24/i),' days','          ', i,' hours'), size=12)
        elif j == 1:
            axs[0,1].plot(p[:,0]/AU,p[:,1]/AU)
            axs[0,1].set_title('N = {}{} {}  $\Delta$t = {}{}' .format(int(365*24/i),' days','          ', i,' hours'), size=12)
        elif j == 2:
            axs[1,0].plot(p[:,0]/AU,p[:,1]/AU)
            axs[1,0].set_title('N = {}{} {}  $\Delta$t = {}{}' .format(int(365*24/i),' days','          ', i,' hours'), size=12)
        elif j == 3:
            axs[1,1].plot(p[:,0]/AU,p[:,1]/AU)
            axs[1,1].set_title('N = {}{} {}  $\Delta$t = {}{}' .format(int(365*24/i),' days','          ', i,' hours'), size=12)
        j += 1
    for ax in axs.flat:
        ax.set(xlabel='x position [AU]', ylabel='y position [AU]')
        ax.label_outer()
    fig.suptitle('Stability of our solution given different timesteps', fontsize=18, y = 0.93)    
    plt.show()
    
    earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,30e3)), 6e24) #re-initializing
    sun = CelestialBody('Sun', np.array((0,0)), np.array((0,0)), 2e30)
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
    plt.title('Total Energy of %s in a complete orbit' %(earth.name), pad=20)
    plt.plot(range(len(ETOT)),ETOT, label='Total Energy')
    plt.xlabel('Days')
    plt.ylabel('Energy [J]')
    plt.legend()
    plt.show()

"""
Solves task 3d
Plots the angular momentum of Earth over one year
Call the function from the console
"""
def proj_3d():
    AU = 149.6e9
    earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,30e3)), 6e24) #re-initializing
    sun = CelestialBody('Sun', np.array((0,0)), np.array((0,0)), 2e30)
    l = earth.AngularMom(sun,365,24*3600)
    plt.figure()
    plt.plot(range(len(l)),l)
    plt.xlabel('Time [Days]')
    plt.ylabel('Angular Momentum [kg m^2/s]')
    plt.title('Angular Momentum for a complete orbit', pad=20)
    plt.show()

"""
Solves task 3e
Firstly plots orbits given different values of beta over 50 years
Secondly plots potential, kinetic and total energy of a Sun-Earth system given 
an elliptical orbit over one year
Thirdly plots the angular momentum of an elliptical orbit over one year.
Lastly, plots the effects of differing beta values on total energy and 
angular momentum over one year.
Call the function from the console
"""
    
def proj_3e():
    AU = 149.6e9
    betas = np.linspace(2,2.0275,12)
    years = 100
    bOrbits = np.zeros((len(betas),years*365,2))
    sun = CelestialBody('Sun', np.array((0,0)), np.array((0,0)), 2e30)
    j = 0
    for i in betas:
        earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,30e3)), 6e24) #re-initializing for each loop
        earth.Verlet(sun,years*365,24*3600,i)
        bOrbits[j] = earth.p
        j += 1
    bOrbits = bOrbits/AU
    plt.figure(figsize=(7,7))
    for i in range(len(betas)):
        b = betas[i]
        plt.plot(bOrbits[i,:,0]/AU, bOrbits[i,:,1]/AU, label='B = %.4f' %(b))
    plt.legend(loc='upper right')
    plt.title('Earths orbit given different values of Beta over 100 years', pad=20)
    plt.xlabel('x position [AU]')
    plt.ylabel('y position [AU]')
    plt.xlim(-4*1e-10,0.25*1e-10)
    plt.ylim(-0.5*1e-10, 1.5*1e-10)
    plt.show()
    
    earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,23.72e3)), 6e24) #re-initializing
    U, K, ETOT = earth.Energies(sun,365,24*3600)
    plt.figure()
    plt.title('All energies of %s in a complete elliptical orbit' %(earth.name), pad=20)
    plt.plot(range(len(ETOT)),U, label='Potential Energy')
    plt.plot(range(len(ETOT)),K, label='Kinetic Energy')
    plt.plot(range(len(ETOT)),ETOT, label='Total Energy')
    plt.xlabel('Time [Days]')
    plt.ylabel('Energy [J]')
    plt.legend()
    plt.show()
    plt.figure()
    
    plt.title('Total Energy of %s in a complete elliptical orbit' %(earth.name), pad=20)
    plt.plot(range(len(ETOT)),ETOT, label='Total Energy')
    plt.xlabel('Time [Days]')
    plt.ylabel('Energy [J]')
    plt.legend()
    plt.show()
    
    plt.figure()
    
    earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,23.72e3)), 6e24) #re-initializing
    l = earth.AngularMom(sun,365,24*3600)
    plt.plot(range(len(l)),l)
    plt.xlabel('Time [Days]')
    plt.ylabel('Angular Momentum [kg m^2/s]')
    plt.title('Angular Momentum for a complete elliptical orbit', pad=20)
    plt.show()
    plt.figure(figsize=(7,7))
    
    j = 0
    energies = np.zeros((len(betas),365*3))
    for i in betas:
        earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,23.72e3)), 6e24) #re-initializing for each loop
        U, K, energies[j] = earth.Energies(sun,3*365,24*3600,i)
        j += 1
    for i in range(len(betas)):
        b = betas[i]
        plt.plot(range(3*365), energies[i], label='B = %.4f' %(b))
    plt.legend()
    plt.xlabel('Time [Days]')
    plt.ylabel('Energy [J]')
    plt.title('Energies for different values of Beta')
    plt.show()
    j = 0
    plt.figure()
    AngMom = np.zeros((len(betas),365*3))
    for i in betas:
        earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,23.72e3)), 6e24) #re-initializing for each loop
        AngMom[j] = earth.AngularMom(sun,3*365,24*3600,i)
        j += 1
    for i in range(len(betas)):
        b = betas[i]
        plt.plot(range(3*365), AngMom[i], label='B = %.4f' %(b))
    plt.legend()
    plt.xlabel('Time [Days]')
    plt.ylabel('Angular Momentum [kg m^2/s]')
    plt.title('Angular Momentum for different values of Beta', pad=20)
    plt.show()
    
"""
Solves task 3f
Calculates the theoretical escape velocity.
After this, creates a plot of Earths Orbit around the sun given intial velocities
from 35000 m/s to 45000 m/s.
Examine whether theoretical is the same as numercial
Call the function from the console
"""
def proj_3f():
    G = 6.67e-11; AU = 149.6e9
    sun = CelestialBody('Sun', np.array((0,0)), np.array((0,0)), 2e30)
    expected_escape_vel = np.sqrt(2*G*sun.mass/AU)
    print('Expected Escape Velocity:', expected_escape_vel)
    vels = np.linspace(35e3,45e3,11)
    plt.figure(figsize=(7.2,7.2))
    for i in vels:
        earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,i)), 6e24) #re-initializing for each loop
        earth.Verlet(sun,75*365,24*3600) #orbit over 75 years
        plt.plot(earth.p[:,0]/AU,earth.p[:,1]/AU, label='Orbit for intial velocity %.0f m/s' %(i))
    plt.legend(loc='upper right')
    plt.xlabel('x position [AU]')
    plt.ylabel('y position [AU]')
    plt.title('Earths orbit around the sun for different initial velocities over 75 years')
    plt.show()



"""
Solves task 3g
Plots three different plots of our three body system with the Earth, the Sun and Jupiter.
The first plot is using real values, the second Jupiter has 10x its real mass 
and in the third Jupiter has 1000x its real mass
Call the function from the console
"""
def proj_3g(): 
    AU = 149.6e9
    earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,30e3)), 6e24) #re-initializing
    sun = CelestialBody('Sun', np.array((0,0)), np.array((0,0)), 2e30) #re-initializing
    jupiter = CelestialBody('Jupiter', np.array((0,5.2*AU)), np.array((13e3,0)), 1.9e27)
    bodies= [earth, jupiter, sun]
    CelestialBody.VerletMultiBody(bodies, 30*365, 24*3600)
    plt.figure()
    plt.plot(earth.p[:,0]/AU,earth.p[:,1]/AU, label='Earth')
    plt.plot(jupiter.p[:,0]/AU, jupiter.p[:,1]/AU, label='Jupiter')
    plt.plot(sun.p[:,0]/AU, sun.p[:,1]/AU, label='Sun')
    plt.xlabel('x position [AU]')
    plt.ylabel('y position [AU]')
    plt.title('Earth and Jupiters Orbit with real Jupiter mass over 30 years')
    plt.legend()
    plt.show()
    
    plt.figure()
    earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,30e3)), 6e24) #re-initializing
    sun = CelestialBody('Sun', np.array((0,0)), np.array((0,0)), 2e30) #re-initializing
    jupiter = CelestialBody('Jupiter', np.array((0,5.2*AU)), np.array((13e3,0)), 10*1.9e27)
    bodies= [earth, jupiter, sun]
    CelestialBody.VerletMultiBody(bodies, 30*365, 24*3600)
    plt.plot(earth.p[:,0]/AU,earth.p[:,1]/AU, label='Earth')
    plt.plot(jupiter.p[:,0]/AU, jupiter.p[:,1]/AU, label='Jupiter')
    plt.plot(sun.p[:,0]/AU, sun.p[:,1]/AU, label='Sun')
    plt.xlabel('x position [AU]')
    plt.ylabel('y position [AU]')
    plt.title('Earth and Jupiters Orbit with 10 x Jupiter mass over 30 years')
    plt.legend()
    plt.show()
    
    plt.figure()
    earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,30e3)), 6e24) #re-initializing
    sun = CelestialBody('Sun', np.array((0,0)), np.array((0,0)), 2e30) #re-initializing
    jupiter = CelestialBody('Jupiter', np.array((0,5.2*AU)), np.array((13e3,0)), 1000*1.9e27)
    bodies= [earth, jupiter, sun]
    CelestialBody.VerletMultiBody(bodies, 30*365, 24*3600)
    plt.plot(earth.p[:,0]/AU,earth.p[:,1]/AU, label='Earth')
    plt.plot(jupiter.p[:,0]/AU, jupiter.p[:,1]/AU, label='Jupiter')
    plt.plot(sun.p[:,0]/AU, sun.p[:,1]/AU, label='Sun')
    plt.xlabel('x position [AU]')
    plt.ylabel('y position [AU]')
    plt.title('Earth and Jupiters Orbit with 1000 x Jupiter mass over 30 years')
    plt.legend()
    plt.show()

"""
Solves task 3h
Creates two plots, the first one is a fifty year plot of our three-body system
where the Sun has had it's initial conditions defined so as to make the total
momentum of the system zero. Additionally the sun has been moved so that the
center of mass in the system is (0,0). Compare the stability of this system to
the "uncorrected" system we create and plot in 3g
The second plot is a full plot of our solar system with all 8 planets over a
time period of 200 years.
Pluto is commented but can be uncommented and included in the bodies list to plot
all 8 planets + Pluto, but be warned that Pluto's orbit time significantly deviated
from it's theoretical orbit time, and was therefore not included.
Call the function from the console
"""
def proj_3h():
    AU = 149.6e9
    comx = (1.9e27*5.2*AU)/(2e30+6e24+1.9e27)
    comy = (6e24*AU)/(2e30+6e24+1.9e27)
    sun = CelestialBody('Sun', np.array((-comx,-comy)), np.array((-12.35,-0.09)), 2e30) #re-initializing
    earth = CelestialBody('Earth', np.array((1*AU,0)), np.array((0,30e3)), 6e24) #re-initializing
    jupiter = CelestialBody('Jupiter', np.array((0,5.2*AU)), np.array((13e3,0)), 1.9e27)
    bodies = [earth,jupiter, sun]
    CelestialBody.VerletMultiBody(bodies, 50*365, 24*3600)
    plt.figure()
    plt.plot(earth.p[:,0]/AU,earth.p[:,1]/AU, label='Earth')
    plt.plot(jupiter.p[:,0]/AU, jupiter.p[:,1]/AU, label='Jupiter')
    plt.plot(sun.p[:,0]/AU, sun.p[:,1]/AU, label='Sun')
    plt.xlabel('x position [AU]')
    plt.ylabel('y position [AU]')
    plt.title('Earth and Jupiters Orbit with real Jupiter mass over 50 years\nwith initialized sun with regards to Center of Mass')
    plt.legend()
    plt.show()
    
    """ center of mass for system with all 8 planets (and Pluto) (antiquated, but 
    kept around as writing out all of this takes a while, and I'd rather not do it
    again should I ever need it)
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
    bodies = [sun,mercury,venus,earth,mars,jupiter,saturn,uranus,neptune] #uncomment above and add pluto in this list to include it in the plot
    CelestialBody.VerletMultiBody(bodies, 200*365, 24*3600)
    plt.figure(figsize=(10,10))
    for body in bodies:
        plt.plot(body.p[:,0]/AU, body.p[:,1]/AU, label='%s' %(body.name))
    plt.xlabel('x position [AU]')
    plt.ylabel('y position [AU]')
    plt.title('Solar System Orbits over 200 years')
    plt.legend(loc='upper left')
    plt.show()


"""
Solves task 3i
First let's our Mercury-Sun system run for 100 years with relativistic force.
Next, let's Mercury keep orbiting until the next time it's at the Perihelion.
It then returns these coordinates which we can use to see the Perihelion has 
moved the theoretical 43 arcseconds over a century.
Call the function from the console
""" 
def proj_3i():
    AU = 149.6e9
    sun = CelestialBody('Sun', np.array((0,0)), np.array((0,0)), 1.989e30) #re-initializing
    mercury = CelestialBody('Mercury', AU*np.array((0.3075,0)), np.array((0,58970.74164)), 3.285e23) #re-initializing
    
    mercury.Verlet(sun,100*365,24*3600, rel=True)
    ps = mercury.p
    radchnge = np.arctan2(ps[87,1], ps[87,0]) #find change per year
    mercury = CelestialBody('Mercury', np.array((mercury.p[100*365-1,0],mercury.p[100*365-1,1])),
                            np.array((mercury.v[100*365-1,0],mercury.v[100*365-1,1])), 3.285e23) #re-initializing from where we left off
    pos, hours = mercury.VerletPerihelion(sun,3600, 1e-20)
    pos = pos/AU
    radtot = (36500+(hours)/24)/88 * radchnge #find change over a century + extra days
    while abs(radtot) > 2*np.pi: #reign in rad to be in [-2pi, 2pi]
        if radtot > 0:
            radtot -= 2*np.pi
        elif radtot < 0:
            radtot += 2*np.pi
            
    lst = []
    for i in range(len(mercury.p)): #getting rid of all zeros in our array
        if np.linalg.norm(mercury.p[i]) != 0:
            lst.append(mercury.p[i])
    ps2 = np.array(lst)
    
    final = 72 #plot only the final 72 values of ps
    
    plt.figure()
    plt.scatter(0,0,s=15, color='r', label='Sun')
    plt.plot(ps[len(ps)-final:,0]/AU,ps[len(ps)-final:,1]/AU, label='Mercury', color='g')
    plt.legend(loc='upper right')
    plt.plot(ps2[:,0]/AU,ps2[:,1]/AU,'-.',color='g', label='Mercury post 100 years')
    plt.plot([ps2[len(ps2)-1,0]/AU,0],[ps2[len(ps2)-1,1]/AU,0], ':', color='r', label='r')
    plt.legend(loc='lower right')
    plt.xlabel('x position [AU]')
    plt.ylabel('y position [AU]')
    plt.title('Mercurys orbit over %i days' %(final+hours/24))
    plt.show()
    
    
    theta = np.arctan2(pos[1],pos[0]) #current angle
    real_deg = theta + radtot #find how much the current angle has deviated from expected value
    print('We find the angle to have changed by %.3f"' %(real_deg/(np.pi/648000)))
    
"""
Solves the entire project in one function.
Allows skipping of 3a as it is essentially the same as 3b, although not using OOP.
Beware that this creates a lot of plots, I'd personally recommend calling
functions one by one.
Call the function from the console
"""
def proj_3_full():
    yn = input(str('Would you like to skip 3a? [Y/N] \n'))
    if yn.lower() == 'y':
        print('--Starting 3a--')
        proj_3a()
        print('--3a finished--')
    print('--Starting 3b--')
    proj_3b()
    print('--3b finished--')
    print('--Starting 3c--')
    proj_3c()
    print('--3c finished--')
    print('--Starting 3d--')
    proj_3d()
    print('--3d finished--')
    print('--Starting 3e--')
    proj_3e()
    print('--3e finished--')
    print('--Starting 3f--')
    proj_3f()
    print('--3f finished--')
    print('--Starting 3g--')
    proj_3g()
    print('--3g finished--')
    print('--Starting 3h--')
    proj_3h()
    print('--3h finished--')
    print('--Starting 3i--')
    proj_3i()
    print('--3i finished--')