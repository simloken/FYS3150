# -*- coding: utf-8 -*-
"""
Standalone program for comparing Jacobi's Method to the Bisection Method for
finding eigenvalues
"""
from functions import tridiag, jacobi_solver, bisect_eigens
import numpy as np
import matplotlib.pyplot as plt
import time
Ns = [5,10,15,20]
t = np.zeros((2,len(Ns),1))
j = 0
for i in Ns:
    tt = time.perf_counter()
    jacobi_solver(tridiag(2,5,i))
    t[0,j,0] = time.perf_counter() - tt
    tt = time.perf_counter()
    bisect_eigens(2,5,i)
    t[1,j,0] = time.perf_counter() - tt
    j += 1
plt.plot(Ns, t[0,:])
plt.plot(Ns, t[1,:])
plt.legend(['Jacobi','Bisection'])
plt.xlabel('Square Matrix Length')
plt.ylabel('Time to calculate eigenvalues [s]') 
plt.title('Time to calculate eigenvalues given a matrix length N')
plt.show()
