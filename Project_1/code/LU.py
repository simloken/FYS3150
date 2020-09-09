# -*- coding: utf-8 -*-
"""
Importing modules
"""
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.linalg as sc
"""
Define Functions
"""
def func(x): #our function
    return(100*np.exp(-10*x))
"""
Define variables
"""
ns = [1,2,3,4] #accepted powers of 10
p = int(input("Please enter the power of 10 from [1,4]: "))
if p not in ns:
    raise ValueError("Not an integer in the list [1,4]")
t0 = time.perf_counter() #start time, the first calculation starts here
n = 10**p
h = 1/(n+1)
hsq = h*h #save FLOPs
k = 0
"""
Define empty arrays
"""
A = np.zeros((n,n))
A[0][0] = 2; A[0][1] = -1
y = np.linspace(0,1,n)
"""
Filling arrays/preparing to solve
"""
for i in range(1,n):
    for j in range(3):
        if j == 0:
            A[i][j+k] = -1
        elif j == 1:
            A[i][j+k] = 2
        elif j ==2:
            if i != n-1:
                A[i][j+k] = -1
    k += 1
    
LU = sc.lu_factor(A)
x = sc.lu_solve((LU), func(y))
t1 = time.perf_counter() - t0 #stop time, find time elapsed
print("Time elapsed: %gs for n = %0.0f" %(t1,n))
file = open('LUData.txt', 'a')
file.write("%0.0f, %g \n" %(n,t1))
file.close()
