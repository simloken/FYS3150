# -*- coding: utf-8 -*-
import numpy as np
from sympy import Matrix
from sympy.utilities import lambdify
from sympy.abc import x
import matplotlib.pyplot as plt
"""
All functions are stored here.
"""
"""
Creates a tridiagonal matrix of size n x n with d along the diagonal and a 
running in parallel along it's sides.
Can handle vector and constant input for d. Not for a (though this could
easily be fixed)
Returns said matrix
"""
def tridiag(a,d,n):
    if hasattr(d, "__len__") == True: #checks if d is an array or not
        A = np.zeros((n,n)) #n x n
        A[0][0] = d[0]; A[0][1] = a
        k = 0
        for i in range(1,n):
            for j in range(3):
                if j == 0:
                    A[i][j+k] = a
                elif j == 1:
                    A[i][j+k] = d[i]
                elif j ==2:
                    if i != n-1:
                        A[i][j+k] = a
            k += 1
        return(A)
    else:
        A = np.zeros((n,n)) #n x n
        A[0][0] = d; A[0][1] = a
        k = 0
        for i in range(1,n):
            for j in range(3):
                if j == 0:
                    A[i][j+k] = a
                elif j == 1:
                    A[i][j+k] = d
                elif j ==2:
                    if i != n-1:
                        A[i][j+k] = a
            k += 1
        return(A)
    
"""
Takes the arguments size, diagonal and it's parallels to find the respective
tridiagonal matrix' eigenvalues
Returns eigenvalues
"""    
def eigenvalues(a,d,n):
    cosine = np.pi/(n+1) #given by project text
    lda = np.zeros(n)
    for i in range(0,n):
        lda[i] = d + 2*a*np.cos(cosine*(i+1)) #given by project text
        
    return lda
"""
Takes a matrix and finds it's biggest squared number Max
Uses the coordinates of Max to find rotation coordinates
Returns x and y coordinates
"""
def max_sqr(mat):
    Max = 0
    N = len(mat)
    for i in range(N):
        for j in range(i+1,N): #check only upper part of diagonal (because of symmetry) to save FLOPs
            sqr = mat[i,j]**2 #save FLOPs
            if  sqr >= Max:
                Max = sqr
                i_Max = i; j_Max = j
    return i_Max, j_Max
"""
Takes a matrix, a vector and two "kick-off" coordinates to rotate a vector vec
given a matrix mat.
Returns a transformed matrix and vector
"""
def rotate(mat,vec,k,l):
    n = len(mat)
    tau = (mat[l,l]-mat[k,k])/(2*mat[k,l]) #given by project text
    if tau >= 0:
       t = 1/(tau+np.sqrt(1+tau**2)) #given by project text
    else:
        t = -1/(-tau+np.sqrt(1+tau**2)) #given by project text
    c = 1/np.sqrt(1+t**2)
    s = t*c
    kmat = mat[k,k]
    lmat = mat[l,l]
    cc = c*c; cs = c*s; ss = s*s #save FLOPs
    mat[k,k] = cc*kmat-2*cs*mat[k,l]+ss*lmat
    mat[l,l] = ss*kmat+2*cs*mat[k,l]+cc*lmat
    mat[k,l] = 0; mat[l,k] = 0
    """
    Rotating and assigning new values
    """
    for i in range(n):
        kvec = vec[i,k]
        lvec = vec[i,l]
        vec[i,k] = c*kvec - s*lvec
        vec[i,l] = c*lvec + s*kvec
        if i not in [k, l]:
            kmat = mat[i,k]; lmat = mat[i,l]
            mat[i,k] = c*kmat - s*lmat
            mat[i,l] = c*lmat + s*kmat
            mat[k,i] = mat[i,k]
            mat[l,i] = mat[i,l]
    return (mat, vec)

"""
Solves our Jacobi Rotation algorithm and returns it's eigenvalues and
eigenvectors
"""
def jacobi_solver(mat):
    tol = 1e-10
    n = len(mat)
    eVec = np.identity(n) #empty identity matrix for filling in eigenvectors
    runs = 0 #tracks runs, used for runs per N measurements
    x,y = max_sqr(mat)
    while mat[x,y]**2 >= tol:
        runs += 1
        mat, eVec = rotate(mat,eVec,x,y)
        x,y = max_sqr(mat)
    eVal = np.zeros(n) #empty matrix for filling in eigenvalues and form D
    for i in range(n):
        eVal[i] = mat[i,i]
    return eVal, eVec, runs
"""
Compares the eigenvectors of the lowqest eigenvalues for a matrix
d a 0 0 0 ...
a d a 0 0 ...
0 a d a 0 ...
.........
of length n
"""
def compare_eigvecs(a,d,n):
    plt.plot(np.sort(np.linalg.eig(tridiag(a,d,n))[1][np.argmin(np.linalg.eig(tridiag(a,d,n))[0])]))
    imin = np.argmin(jacobi_solver(tridiag(a,d,n))[0])
    plt.plot(range(len(jacobi_solver(tridiag(a,d,n))[0])),np.sort(jacobi_solver(tridiag(a,d,n))[1][imin]))
    plt.legend(['Analytical','Jacobi'])
    plt.xlabel('Eigenvector number')
    plt.ylabel('Eigenvector length')
    plt.title('Eigenvector comparisons between the Analytical and the Jacobian solutions')
    plt.show()
"""
Takes a square matrix length n and some random constant c.
c represents rho_N - rho_0 given in the assignment.
Different values of c give different accuracies of lambda, should be adjusted
accordingly (trial and error).
"""
def quantumdots_one(n,c): 
    h = c/n
    lgth = np.array(range(1,n+1)) #length array
    rho = lgth*h #rhos
    pot = rho**2 #the potential
    d = 2/(h**2) + pot #given by project text
    e = -1/(h**2) #given by project text
    A = tridiag(e,d,n)
    eVal, eVec, runs = jacobi_solver(A) #retrieve eigenvalues/vectors
    couple = eVal.argsort() #couple and sort them accordingly
    eVal = eVal[couple]
    eVec = eVec[:, couple] 
    for i in range(0,3):
        plt.plot(rho, eVec[:, i])
        print('Lambda%i = %0.4f' %(i+1, eVal[i]))
    plt.xlabel(r'$\rho$')
    plt.ylabel('Solution')
    plt.title('Radial Solution given a constant %i for the first 3 eigenstates \n with N = %i integration points' %(c, n))
    plt.legend(['$\lambda_1$','$\lambda_2$','$\lambda_3$'])
    plt.show()
    
"""
Performs a bisection, looks for roots between a and b for a function f over N
iterations
"""
def bisection(f,a,b,N):
    for n in range(1,N+1):
        c = (a + b)/2
        fc = f(c)
        if f(a)*fc < 0:
            a = a
            b = c
        elif f(b)*fc < 0:
            a = c
            b = b
        elif fc == 0:
            return c
        else:
            return None
    return (a + b)/2   
"""
Finds the characteristic polynomial and parses it as readable lambda function
of variable x
"""
def poly(mat):               
    A = Matrix(mat)
    A = A.charpoly(x='x').as_expr()
    A = lambdify(x,A)
    return A
"""
Finds the eigenvalues given a matrix with d diagonal, a subdiagonal and n length
Defaults to checking for eigenvalues in [0,20].
Not compatible with negative eigenvalues unless neg = True.
Slower than Jacobi :(
"""            
def bisect_eigens(a,d,n,u=10, neg=False):
    L = []
    l = 0
    f = (poly(tridiag(a,d,n)))
    if neg == True:
        l = -u
    for i in np.linspace(l,2*u,10*u):
        eig = bisection(f,i,2*i,500)
        if eig != None:
            eig = round(eig,8) #rounds to eight decimals so as to not get duplicates
            if eig not in L and np.allclose(L,eig,rtol = 1e-3, atol=1e-2) == False: #eliminates duplicates
                L.append(eig)
            elif not L: #first time
                L = [eig]
    L.sort()
    return(L)