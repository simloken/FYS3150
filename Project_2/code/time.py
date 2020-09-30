# -*- coding: utf-8 -*-
"""
Standalone program for calculating runtime for three different algorithms for
finding eigenvalues
"""
from functions import tridiag, eigenvalues, jacobi_solver
import numpy as np
import matplotlib.pyplot as plt
import time
Ns = [5,10,30,50]
t = np.zeros((3,len(Ns),1))
j = 0
for i in Ns:
    tt = time.perf_counter()
    eigenvalues(-5,10,i)
    t[0,j,0] = time.perf_counter() - tt
    tt = time.perf_counter()
    np.linalg.eig(tridiag(-5,10,i))
    t[1,j,0] = time.perf_counter() - tt
    tt = time.perf_counter()
    jacobi_solver(tridiag(-5,10,i))
    t[2,j,0] = time.perf_counter() - tt
    j += 1
plt.plot(Ns, t[0,:])
plt.plot(Ns, t[1,:])
plt.plot(Ns, t[2,:])
plt.legend(['Analytical','Numpy','Jacobi'])
plt.xlabel('Square Matrix Length')
plt.ylabel('Time to calculate eigenvalues [s]') 
plt.title('Time to calculate eigenvalues given a matrix length N')
plt.show()