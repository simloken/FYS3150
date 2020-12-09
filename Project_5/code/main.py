import numpy as np
import matplotlib.pyplot as plt
from functions import G
from classes import Simulation

#np.random.seed(98) # good seed for MCS = 100000
#np.random.seed(33) # good seed for MCS = 10000
sim = Simulation(500, 10000, 7000000, plotting=True)
A,B,C,M = sim.run()
   
M = np.asarray(M)
plt.figure()
plt.plot(range(len(M)), M)
plt.xlabel('Monte Carlo Cycles')
plt.ylabel('Magnetization m')
plt.ylim([-1.05,1.05])
plt.show()
dtlen = int(len(M)/2)
dts = np.linspace(1, dtlen, dtlen)
g = G(dts, M)
plt.figure()
plt.plot(dts, g)
plt.title('Auto correlation')
plt.xlabel('$\Delta t$')
plt.ylabel('$G(\Delta t)$')
plt.show()


def decisionTimes(runs, N, NN, MCS, plotting=False): #pass an object instead? or dont
    decisionTime = []
    times = []
    adt = []
    for i in range(runs):
        
        sim = Simulation(N, NN, MCS, plotting)
        sim.run()
        
        #for j in range(sim.time):
            #decisionTime[j] = (sim.dT[j])
    
        decisionTime.append(sim.dT)
        adt.append(sim.adT)
    lst = []
    for i in decisionTime:
        for j in i:
            lst.append(j)
    lst2 = []        
    for i in adt:
        for j in i:
            lst2.append(j)
        
    dt = np.asarray(lst).ravel()
    adt = np.asarray(lst2).ravel()
    
    return dt, adt

dt, adt = decisionTimes(1, 500, 10000, 7000000)
x = np.linspace(2,10000, 9998)
def powerlaw(x, p, c=1):
    return c*x**(p)
plt.figure()
power = -4/3

plt.loglog(range(len(dt)), 0.1/dt, '+', color='g')
plt.loglog(x, powerlaw(x, (power)), color='r')
plt.title('Correlation between a power law with exponent ~ %g and decision times' %(power))
plt.xlabel('$\\tau$')
plt.ylabel('$P(\\tau)$')
plt.show()
"""
plt.figure()
plt.loglog(range(len(adt)), 0.4/adt, '+', color='g')
plt.loglog(x, powerlaw(x, (power)), color='r')
plt.title('Correlation between a power law with exponent ~ %g and decision times' %(power))
plt.xlabel('$\\tau$')
plt.ylabel('$P(\\tau)$')
plt.show()
"""

def clustered(cB, N, NN, MCS, random=False, plotting=False):
    
    sim = Simulation(N, NN, MCS, cB, clusterd=True, random=random, plotting=plotting)
    A,B,C,M = sim.run()
    plt.figure()
    plt.plot(range(len(M)), M)
    plt.show()
    
    dt = sim.dT
    dt = np.asarray(dt)
    adt = sim.adT
    adt = np.asarray(adt)
    
    return dt, adt
    
dt2, _ = clustered(0.2, 500, 10000, 7000000, random=True, plotting=False)
plt.loglog(x, powerlaw(x, (-3/2)), color='r')
plt.loglog(range(len(dt2)), 0.1/dt2, '+', color='g')
plt.show()


def probas(N, NN, MCS, switch_proba):
    cB=0
    sim = Simulation(N, NN, MCS, cB=cB, clusterd=True, random=False, plotting=False, switch_proba=switch_proba)
    A,B,C,M = sim.run()
    plt.figure()
    plt.title('Magnetization with a noise p=%g' %(switch_proba))
    plt.plot(range(len(M)), M)
    plt.ylim([-1.05,1.05])
    plt.show()
    

probas(500, 10000, 7000000, 3e-6)
    
    