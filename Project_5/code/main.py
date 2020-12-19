import numpy as np
import matplotlib.pyplot as plt
from functions import G, powerlaw
from matplotlib.ticker import ScalarFormatter
from classes import Simulation

"""
    Solves part a), b) and c) of the project.
    N : integer
    The number of elements in the chain
    NN : integer
    Number of places to sample to find spread
    MCS : integer
    Number of Montecarlo cycles.
    plotting : boolean
    Whether or not to plot the opinion spread underway
    
    Returns a bunch of plots, including progress
    
"""

def partABC(N, NN, MCS, plotting=True):
    sim = Simulation(N, NN, MCS, plotting=plotting)
    A,B,C,M = sim.run()
       
    M = np.asarray(M)
    plt.figure()
    plt.plot(range(len(M)), M)
    plt.xlabel('Monte Carlo Cycles')
    plt.ylabel('Magnetization m')
    plt.title('The magnetization of Monte Carlo cycles')
    plt.ylim([-1.05,1.05])
    plt.ticklabel_format(axis='x', style='sci', scilimits=[0,0])
    plt.show()
    dtlen = int(len(M)/2)
    dts = np.linspace(1, dtlen, dtlen)
    g = G(dts, M)
    plt.figure()
    plt.plot(dts, g)
    plt.title('Auto correlation')
    plt.xlabel('$\Delta t$')
    plt.ylabel('$G(\Delta t)$')
    plt.ticklabel_format(axis='x', style='sci', scilimits=[0,0])
    plt.show()




"""
    Studies the correlation between decision times and a power law logarithmicly
    N : integer
    The number of elements in the chain
    NN : integer
    Number of places to sample to find spread
    MCS : integer
    Number of Montecarlo cycles.
    plotting : boolean
    Whether or not to plot the opinion spread underway
    
    Returns decision times in an array
    
"""

def decisionTimes(N, NN, MCS, plotting=False): #pass an object instead? or dont
    decisionTime = []
    adt = []

    
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


"""
    Studies the effects of different concentrations of B
    cB : float
    Percentage of population to have opinion B. Must be within [0, 1] 
    N : integer
    The number of elements in the chain
    NN : integer
    Number of places to sample to find spread
    MCS : integer
    Number of Montecarlo cycles.
    random : boolean
    Whether or not to randomly spread the B opinions or cluster them.
    plotting : boolean
    Whether or not to plot the opinion spread underway
    
    Returns plots of m and decision times
    
"""

def clustered(cB, N, NN, MCS, random=False, plotting=False):
    
    sim = Simulation(N, NN, MCS, cB, clusterd=True, random=random, plotting=plotting)
    A,B,C,M = sim.run()
    plt.figure()
    plt.plot(range(len(M)), M)
    plt.xlabel('Monte Carlo Cycles')
    plt.ylabel('Magnetization m')
    plt.ylim([-1.05,1.05])
    plt.ticklabel_format(axis='x', style='sci', scilimits=[0,0])
    plt.show()
    
    dt = sim.dT
    dt = np.asarray(dt)
    adt = sim.adT
    adt = np.asarray(adt)
    
    return dt, adt
    



"""
    Studies the effects of a noise probability
    N : integer
    The number of elements in the chain
    NN : integer
    Number of places to sample to find spread
    MCS : integer
    Number of Montecarlo cycles.
    switch_proba : None or float
    If not None, then this is the probability that a random element i will flip at a given timestep.
    
    Returns a plot of m
    
"""

def probas(N, NN, MCS, switch_proba):
    cB=0
    sim = Simulation(N, NN, MCS, cB=cB, clusterd=True, random=False, plotting=False, switch_proba=switch_proba)
    A,B,C,M = sim.run()
    plt.figure()
    plt.title('Magnetization with a noise p=%g' %(switch_proba))
    plt.xlabel('Monte Carlo Cycles')
    plt.ylabel('Magnetization m')
    plt.plot(range(len(M)), M)
    plt.ylim([-1.05,1.05])
    plt.ticklabel_format(axis='x', style='sci', scilimits=[0,0])
    plt.show()

def decisionProbas():
    cB = 0
    switch_probas = [2e-4, 2e-3, 2e-2, 2e-1]
    for i in switch_probas:
        sim = Simulation(500, 10000, 2000000, cB=cB, clusterd=True, random=False, plotting=False, switch_proba=i)
        sim.run()
        adt = sim.adT
        adt = np.asarray(adt)
        plt.loglog(range(len(adt)), 3/adt, '+', label='p: %g ' %(i))
        
    plt.xlabel('$\\tau$')
    plt.ylabel('$P(\\tau)$')
    plt.title('Decision times for different noises p')
    plt.legend()
    plt.show()


#part ABC    
partABC(500, 10000, 10000000)

#decision times
dt, adt = decisionTimes(500, 10000, 10000000)
x = np.linspace(0,100000, 100001)

plt.figure()
power = -3/2

plt.loglog(range(len(adt)), 3/adt, '+', color='g') #3 here is just a scaling parameter, else the plot goes haywire.
plt.loglog(x, powerlaw(x, (power)), color='r')
plt.title('Correlation between a power law with exponent ~ %g and decision times' %(power))
plt.xlabel('$\\tau$')
plt.ylabel('$P(\\tau)$')
plt.show()

#clustered B's
dt2, adt2 = clustered(0.2, 500, 10000, 10000000, random=True, plotting=False)
plt.loglog(x, powerlaw(x, (-3/2)), color='r')
plt.loglog(range(len(adt2)), 3/adt2, '+', color='g')
plt.title('Correlation between a power law with exponent ~ %g and decision times' %(power))
plt.xlabel('$\\tau$')
plt.ylabel('$P(\\tau)$')
plt.show()

#probability and noise
probas(500, 10000, 10000000, 8e-6)
    
#decisiontimes for different probabilities
decisionProbas()
