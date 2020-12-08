import numpy as np

def G(dt, m): #this doesnt work
    g = np.zeros(len(dt))
    mm = np.mean(m)
    M = np.sum(m - mm)
    up = M
    low = (M)**2
    for i in range(len(g)):
        g[i] = (m[i]-mm)*up/low
        
    return g