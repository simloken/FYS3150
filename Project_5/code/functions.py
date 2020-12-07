import numpy as np

def G(m): #this doesnt work
    g = np.zeros(len(m)-1)
    for i in range(len(m)-1):
        mm = np.mean(m)
        up = np.sum((m[i]-mm)*(m[i+1]-mm))
        low = (np.sum(m[i]-mm))**2
        g[i] = up/low
    return g