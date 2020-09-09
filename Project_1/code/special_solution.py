# -*- coding: utf-8 -*-
"""
Importing modules
"""
import numpy as np
import matplotlib.pyplot as plt
import time
"""
Define Functions
"""
def aU(x): #analytical solution given by the exercise
    return(1-(1-np.exp(-10))*x-np.exp(-10*x))
def func(x): #our function
    return(100*np.exp(-10*x))
def spec(u,f,n,d_t,f_t):
    f_t[1] = f[1]
    for i in range(2,n+1): #forward sub
        f_t[i] = f[i]+(f_t[i-1])/d_t[i-1]
    u[n] = (f_t[n])/d_t[n]
    for j in range(n-1,0,-1): #backwards sub
        u[j]= (f_t[j]+u[j+1])/d_t[j]
    return u

"""
Define variables
"""
ns = [1,2,3,4,5,6,7] #accepted powers of 10
p = int(input("Please enter the power of 10 from [1,7]: "))
if p not in ns:
    raise ValueError("Not an integer in the list [1,7]")
t0 = time.perf_counter() #start time, the first calculation starts here
n = 10**p
h = 1/(n+1)
hsq = h*h #save FLOPs

"""
Define empty arrays
"""
x = np.zeros(n+2); f = np.zeros(n+2); a = np.zeros(n+2) #arrays used to solve
d_t = np.zeros(n+1); f_t = np.zeros(n+1); #tilde arrays
u = np.zeros(n+2); y = np.linspace(0,1,n+2) #num.sol. array and filler array

"""
Filling arrays/preparing to solve
"""
for i in range(n+2):
    x[i] = h*i
    f[i] = func(x[i])*hsq
    a[i] = aU((x[i]))
    
for j in range(1,n+1):
    d_t[j] = (j+1)/j
"""
Solving
"""
rU = spec(u,f,n,d_t,f_t) #numerical solution is calculated here
plt.plot(y,aU(y), '-b') #plot analytical
plt.plot(y,rU, '--r') #plot numerical
plt.legend(['Analytical Solution', 'Numerical Solution'])
plt.grid(b=True) #grids
t1 = time.perf_counter() - t0 #stop time, find time elapsed
plt.show()
print("Time elapsed: %gs for n = %0.0f" %(t1,n))
