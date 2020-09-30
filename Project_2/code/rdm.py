# -*- coding: utf-8 -*-
from functions import tridiag, eigenvalues, max_sqr, rotate, jacobi_solver, bisection, poly, bisect_eigens, compare_eigvecs, quantumdots_one
import time
import matplotlib.pyplot as plt
import numpy as np
#eigenvalues printed
"""
print('Jacobi ',np.sort(jacobi_solver(tridiag(2,5,5))[0]))
print('Analytical ',np.sort(eigenvalues(2,5,5)))
print('Numpy ', np.sort(np.linalg.eig(tridiag(2,5,5))[0]))
"""
#time taken to bisect
"""
t = time.perf_counter()
B = bisect_eigens(2,5,5)
t1 = time.perf_counter() - t
print(B)
print(t1)
"""
#plotting eigenvectors for all algorithms
"""
n = 10
u = np.zeros(n)
for i in range(n):
    u[i] = np.sin((i*np.pi)/n)
v = jacobi_solver(tridiag(2,5,n))[1][np.argmin(jacobi_solver(tridiag(2,5,n))[0])]
w = np.linalg.eig(tridiag(2,5,n))[1][np.argmin(np.linalg.eig(tridiag(2,5,n))[0])]
plt.plot(range(len(u)),u)
plt.plot(range(len(v)),v)
plt.plot(range(len(w)),w)
plt.show()
"""
#plotting eigenvalues for all algorithms
"""
n = 100
u = np.sort(jacobi_solver(tridiag(2,5,n))[0])
v = np.sort(eigenvalues(2,5,n))
w = np.sort(np.linalg.eig(tridiag(2,5,n))[0])
plt.plot(range(len(u)),u)
plt.plot(range(len(v)),v)
plt.plot(range(len(w)),w)
plt.show()
"""
#2d
"""
quantumdots_one(200,10)
"""