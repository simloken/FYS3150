import numpy as np
import matplotlib.pyplot as plt
from functions import G
from classes import Simulation

#np.random.seed(98) # good seed for MCS = 100000
#np.random.seed(33) # good seed for MCS = 10000
sim = Simulation(1000, 10000, 10000, plotting=True)
A,B,C,M = sim.run()
plt.figure()
plt.plot(range(len(M)), M)
plt.show()
#g = G(M)
#plt.figure()
#plt.plot(range(len(g)), g)
#plt.show()


def decisionTimes(runs, N, NN, MCS, plotting=False):
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

#dt = decisionTimes(1, 1000, 10000, 100000)

#plt.loglog(np.arange(1, 2000, 1), np.arange(1, 2000, 1) ** (-3/2), color='r')

#plt.plot(dt, dt, '+', color='g')
#plt.show()