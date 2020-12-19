import numpy as np

"""
   The autocorrelation function
   dt : array-like
   Tells the return array how long it must be
   m : array-like
   The magnetization
"""
def G(dt, m): #this doesnt work? This is probably wrong.
    g = np.zeros(len(dt))
    mm = np.mean(m)
    M = np.sum(m) - mm
    up = M
    low = (M)**2
    for i in range(len(g)):
        g[i] = (M*m[i]-mm)*up/low
        
    return g

"""
    The powerlaw function
    x : array-like
    Input array to study over, x-coordinates essentially
    p : float
    The power of x
    c : float
    A constant.
"""
def powerlaw(x, p, c=10):
    return c*x**(p)