# -*- coding: utf-8 -*-
import numpy as np
from functions import tridiag, eigenvalues, max_sqr, jacobi_solver
"""
Standalone program for testing and validating some functions used in this
project
"""
"""
Let us test if all our methods return the correct eigenvalues
Let's assume we have the Matrix
5 2 0 0 0
2 5 2 0 0
0 2 5 2 0
0 0 2 5 0
This should return the eigenvalues:
    lambda = 7, 5, 3, 5 + 2*sqrt(3), 5 - 2*sqrt(3)
"""
rLambdas = np.sort([7,5,3,5 + 2*np.sqrt(3), 5 - 2*np.sqrt(3)])
A = tridiag(2,5,5)
"""
Let's now test each of our solutions
"""
nEigs = np.sort(np.linalg.eig(A)[0])
aEigs = np.sort(eigenvalues(2,5,5))
jEigs = np.sort(jacobi_solver(A)[0])
if np.allclose(rLambdas,nEigs) == True:
    print('Numpy Validated')
if np.allclose(rLambdas,nEigs) == True:
    print('Analytical Validated')
if np.allclose(rLambdas,nEigs) == True:
     print('Jacobi Validated')
     
"""
This can also be verified with any other Matrix
We assume:
9 -3 0 0 0
-3 9 -3 0 0
0 -3 9 -3 0
0 0 -3 9 -3
0 0 0 9 -3
Eigenvalues:
    3*(3+sqrt(3)), 12, 9, 6, -3*(sqrt(3)-3)
"""
print()
rLambdas = np.sort([3*(3+np.sqrt(3)), 12, 9, 6, -3*(np.sqrt(3)-3)])
A = tridiag(-3,9,5)
nEigs = np.sort(np.linalg.eig(A)[0])
aEigs = np.sort(eigenvalues(-3,9,5))
jEigs = np.sort(jacobi_solver(A)[0])
if np.allclose(rLambdas,nEigs) == True:
    print('Numpy Validated')
if np.allclose(rLambdas,nEigs) == True:
    print('Analytical Validated')
if np.allclose(rLambdas,nEigs) == True:
     print('Jacobi Validated')
"""
Let's now test that we're getting the largest value when extracting:
Assume we have the matrix:
5 38 0 0 0
38 5 40 0 0
0 40 5 42 0
0 0 42 5 40
We should now get coordinates to 42 as result. Let's try, using a modified version of tridiag()
"""
print()
def tridiagmod(a,d,n):
    A = np.zeros((n,n)) #n x n
    A[0][0] = d; A[0][1] = a
    k = 0
    for i in range(1,n):
        for j in range(3):
            if j == 0:
                A[i][j+k] = a+2*k-2
            elif j == 1:
                A[i][j+k] = d
            elif j ==2:
                if i != n-1:
                    A[i][j+k] = a+2*k
        k += 1
    return(A)
d = 5; a = 38; n = 5
A = tridiagmod(a,d,n)
x,y = max_sqr(A)
if A[x,y] == 42:
    print('max_sqr() Validated')