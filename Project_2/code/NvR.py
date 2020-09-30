# -*- coding: utf-8 -*-
"""
Standalone program for plotting N versus Rotations
Runtime is pretty long, beware! The "soft-cap" seems to be at around N = 100,
so remove N = 150 if its taking too long
"""
from functions import tridiag, jacobi_solver
import numpy as np
import matplotlib.pyplot as plt
Ns = [10,20,30,50,100,150] #N values
data = np.zeros((len(Ns),2))
j = 0
for i in Ns:
    data[j][0] = i #x-axis
    data[j][1] = jacobi_solver(tridiag(-1,2,i))[2] #y-axis, returns iterations
    j += 1
plt.xlabel('Square Matrix Length')
plt.ylabel('Number of iterations needed')
plt.title('Iterations needed with the Jacobi Rotation Algorithm given a length N')
plt.plot(data[:,0], data[:,1])
plt.show()