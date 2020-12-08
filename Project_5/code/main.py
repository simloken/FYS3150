import numpy as np
import matplotlib.pyplot as plt
from functions import G
from classes import Simulation
from multiclass import Simulation as Simp

#np.random.seed(98) # good seed for MCS = 100000
#np.random.seed(33) # good seed for MCS = 10000
sim = Simulation(1000, 10000, 10000, plotting=True)
#sim = Simp(1000, 10000, 10000, plotting=True)
A,B,C,M = sim.run()
   

plt.figure()
plt.plot(range(len(M)), M)
plt.show()
#g = G(M)
#plt.figure()
#plt.plot(range(len(g)), g)
#plt.show()


def decisionTimes(runs, N, NN, MCS, plotting=False): #pass an object instead?
    decisionTime = []
    times = []
    for i in range(runs):
        
        sim = Simulation(N, NN, MCS, plotting)
        sim.run()
        
        #for j in range(sim.time):
            #decisionTime[j] = (sim.dT[j])
    
        decisionTime.append(sim.dT)
        
    lst = []
    for i in decisionTime:
        for j in i:
            lst.append(j)
            
        
    dt = np.asarray(lst).ravel()
    
    return dt

dt = decisionTimes(1, 1000, 10000, 10000)
x = np.linspace(2,10000, 9998)
def power(x, p):
    return x**(p)
plt.loglog(x, power(x, (-3/2)), color='r')
plt.loglog(range(len(dt)), 0.1/dt, '+', color='g')
plt.show()
weights=np.ones_like(dt)/len(dt)
plt.plot(x, power(x, (-1/3)))
plt.hist(dt, bins=10, weights=weights, histtype='bar', ec='black', color='r')
plt.show()



def clustered(cB, N, NN, MCS, random=False, plotting=False):
    
    sim = Simulation(N, NN, MCS, cB, clusterd=True, random=random, plotting=plotting)
    A,B,C,M = sim.run()
    plt.figure()
    plt.plot(range(len(M)), M)
    plt.show()
    
    dt = sim.dT
    dt = np.asarray(dt)
    
    return dt
    
dt2 = clustered(0.5, 1000, 10000, 10000, random=True, plotting=True)
plt.loglog(x, power(x, (-3/2)), color='r')
plt.loglog(range(len(dt2)), 0.1/dt2, '+', color='g')
plt.show()