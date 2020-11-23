import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
"""

    Periodic boundary function. Called in monteCarlo()
    i : integer
    Either an x or y coordinate
    mx : integer
    The max allowed value in our lattice, that is, L for a lattice L x L.
    add : integer
    Parameter added in order to handle boundary cases.
    
    Used everywhere
"""
def period(i, mx, add):
    return (i+mx+add)%mx

"""
    The bread and butter of this project. Performs a number of Monte Carlo sweeps
    across our 2D ising model, then returns varying values depending on input
    booleans
    
    T : float
    The temperature of our material
    L : integer
    The dimension of our lattice, L x L
    cycles : integer
    The number of Monte Carlo sweeps to do
    progress : boolean
    Tells the function whether or not to print a status update for how far along
    it is in it's calculation every set interval
    state : boolean
    Tells the function whether or not to randomize the initial lattice. If False
    then the lattice is uniform
    accept : boolean
    Tells the function whether or not to also return a counter for how many 
    accepted configurationswe got
    proba : boolean
    Tells the function if we're interested in finding out the probability of
    different energies. If this is True, only the energy array is returned
    
    Used pretty much everywhere, with the booleans modifying its usage.
"""
def monteCarlo(T, L, cycles, progress=True, state=False, accept=False, proba=False):
    if progress==True:
        print('ID %i working through %i cycles' %(int(mp.current_process()._identity[0]), cycles))
    if state == True:
        spintrix = np.random.choice([-1,1], size=(L,L)) #either up or down
    else:
        spintrix = np.full((L,L), fill_value=1) #all same
        
    nonTen = False
    count = 0
    probaE = np.zeros(cycles)
    if cycles % 10 != 0:
        nonTen = True
    E = 0; M = 0
    Eex = 0; E2ex = 0
    Mex = 0; M2ex = 0; absMex = 0
    
    w = np.zeros(17,np.float64)
    for de in range(-8,9,4):
        w[de+8] = np.exp(-de/T)
    
    M = spintrix.sum()
    for j in range(L): 
        for i in range(L):
            E -= spintrix.item(i,j)*\
                 (spintrix.item(period(i,L,-1),j) + spintrix.item(i,period(j,L,1)))       
    k = 0             
    for i in range(cycles):
        for s in range(L**2):
            x = np.random.randint(0,L)
            y = np.random.randint(0,L)
            deltaE = 2*spintrix.item(x,y)*(
                     (spintrix.item(period(x,L,-1), y) +
                      spintrix.item(period(x,L,1),  y) +
                      spintrix.item(x, period(y,L,-1)) +
                      spintrix.item(x, period(y,L,1))))
            if np.random.random() <= w[deltaE+8]:
                spintrix[x,y] *= -1
                M += 2*spintrix[x,y]
                E += deltaE
                count += 1
        probaE[i] += E
        Eex    += E
        E2ex   += E**2
        Mex    += M
        M2ex   += float(M**2)
        absMex += int(np.fabs(M))
       
            
     
        if progress==True:
            if nonTen == True:
                if i % int(cycles/10) == 0 and i != 0: #special case so we don't get percentages higher than 100
                    k += 1
                    print('ID %i at %i%%' %(int(mp.current_process()._identity[0]), int(np.round(k*9.9))))
            else:
                if i % int(cycles/10) == 0:
                    k += 1
                    print('ID %i at %i%%' %(int(mp.current_process()._identity[0]), int(k*10)-10))             
    store = np.array((T, Eex, E2ex, Mex, M2ex, absMex))
    
    if progress==True:
        print('ID %i at %i%%' %(int(mp.current_process()._identity[0]), int(100)))
    if proba == True:
        return probaE
    if accept == True:
        return store, count
    else:
        return store
"""
    Our sorting function. Has two different algorithms, depending on whether
    the boolean alt is True or not.
    
    arr : array-like
    Our input array with numerical values
    T : array-like or float
    The temperature of our material
    L : integer
    Dimension of our lattice, L x L
    cycles : array-like or integer
    Number of Monte Carlo cycles to generate our numerical data.
    extras : boolean
    Dictates whether or not to also calculate heat capacity and magnetic
    susceptability
    alt : boolean
    Tells us which input is a variable. If False, then it's assumed that 
    T is an array. If True, then we assume that cycles is an array-like.
    
    Used everywhere.
"""
def sorter(arr,T, L, cycles, extras=False, alt=False):
    if alt == False:
        sarr = arr[np.argsort(arr[:,0])] #sort by temperature
        TL = len(T)
        narr = np.zeros((TL, 6))
        for i in range(TL):
            narr[i][0] = T[i]
        for k in range(TL):
            for i in sarr:
                if i[0] == narr[k][0]:
                    for j in range(len(i) -1):
                        narr[k][j+1] += i[j+1]
        sarr = np.zeros((TL, 6))
        for i in range(TL):                
            sarr[i] = normer(narr[i], L, cycles)  
        if extras == True:
            Cv = np.zeros(TL); X = np.zeros(TL)
            for i in range(TL):
                Cv[i] = CV(sarr[i,2], T[i])
                X[i] = chi(sarr[i,5], T[i])
            
            return sarr, Cv, X
        return sarr
    elif alt == True:

        if arr.shape == (6,):
            sarr= normer(arr, L, cycles)
        else:
            narr = np.sum(arr, axis=0)
            sarr= normer(narr, L, cycles)
            sarr[0] = T
        return sarr
"""
    A normalizing function. Never called on its own but instead called in sorter()
    
    arr : array-like
    Our input array with numerical values
    L : integer
    Dimension of our lattice, L x L
    cycles : integer
    Number of Monte Carlo cycles to generate our numerical data.
    
    Used everywhere.
"""    
def normer(arr, L, cycles):
    T = arr[0]
    Eex = arr[1]; E2ex = arr[2]; Mex = arr[3]; M2ex = arr[4]; absMex = arr[5]
    
    Eex       /= float(cycles)
    E2ex      /= float(cycles)
    Mex       /= float(cycles)
    M2ex      /= float(cycles)
    absMex    /= float(cycles)
    
    L2 = L*L
    E_variance  = (E2ex-Eex**2)/float(L2*T*T)
    M_variance  = (M2ex-absMex**2)/float(L2*T)
    
    Eex       /= float(L2)
    Mex       /= float(L2)
    absMex    /= float(L2)
    
    narr = np.array((T, Eex, E_variance, Mex, M_variance, absMex))
    return narr

"""
    A function that handles a set of inputs to plot the stabilization/equilibrium
    of our numerical data for differing number of Monte Carlo cycles
    
    arr : array-like
    Input array with our numerical values
    cL : array-like
    Array-like with all different cycle values. Values must be integers
    L : integer
    Dimension of our lattice, L x L
    T : float
    Temperature of our material
    w : string
    Parses into plotting title whether the data was generated using an ordered or
    disordered/random initial lattice.
    ana : boolean
    Whether or not we want to include analytical values in our plot.
    Should only be true for lattices of dimension 2 x 2.
    
    Used in part c) and d)
"""
def Stabilize(arr, cL, L, T, w, ana=True):
    if ana==True:
        analytical = analytics(T)/(L**2)
        plt.figure()
        plt.semilogx(cL, analytical[0]*np.ones(len(cL)),'r.-', label='Analytic Mean E')
        plt.semilogx(cL, arr[:,1], 'gx-', label='Num. Mean E')
        plt.title('Expectation E against analytical result for T=%g \n%s ordering' %(T, w))
        plt.xlabel('Monte Carlo cycles')
        plt.legend()
        plt.show()
        plt.figure()
        plt.semilogx(cL, analytical[1]*np.ones(len(cL)), 'r.-', label='Analytic Mean M')
        plt.semilogx(cL, arr[:,5], 'gx-', label='Num.Mean M')
        plt.title('Expectation M against analytical result for T=%g \n%s ordering' %(T, w))
        plt.xlabel('Monte Carlo cycles')
        plt.legend()
        plt.show()
        plt.figure()
        plt.semilogx(cL, analytical[2]*np.ones(len(cL)), 'r.-', label='Analytic Specific Heat')
        plt.semilogx(cL, CV(arr[:,2], T), 'gx-', label='Num. Specific Heat')
        plt.title('Heat Capacity against analytical result for T=%g \n%s ordering' %(T, w))
        plt.xlabel('Monte Carlo cycles')
        plt.legend()
        plt.show()
        plt.figure()
        plt.semilogx(cL, analytical[3]*np.ones(len(cL)),'r.-', label='Analytic Mag. Susc.')
        plt.semilogx(cL, chi(arr[:,4], T), 'gx-', label='Num. Mag. Susceptibility')
        plt.title('Magnetic Susceptibility against analytical result for T=%g \n%s ordering' %(T, w))
        plt.xlabel('Monte Carlo cycles')
        plt.legend()
        plt.show()
    else:
        plt.figure()
        plt.semilogx(cL, arr[:,1], 'gx-', label='Num. Mean E')
        plt.title('Expectation E for T=%g \n%s ordering' %(T, w))
        plt.xlabel('Monte Carlo cycles')
        plt.legend()
        plt.show()
        plt.figure()
        plt.semilogx(cL, arr[:,5], 'gx-', label='Num.Mean M')
        plt.title('Expectation M for T=%g \n%s ordering' %(T, w))
        plt.xlabel('Monte Carlo cycles')
        plt.legend()
        plt.show()
        plt.figure()
        plt.semilogx(cL, CV(arr[:,2], T), 'gx-', label='Num. Specific Heat')
        plt.title('Heat Capacity for T=%g \n%s ordering' %(T, w))
        plt.xlabel('Monte Carlo cycles')
        plt.legend()
        plt.show()
        plt.figure()
        plt.semilogx(cL, chi(arr[:,4], T), 'gx-', label='Num. Mag. Susceptibility')
        plt.title('Magnetic Susceptibility for T=%g \n%s ordering' %(T, w))
        plt.xlabel('Monte Carlo cycles')
        plt.legend()
        plt.show()

"""
    The analytical value function. Takes a temperature T and returns it's corresponding
    analytical values. Keep in mind that it presupposes that the coupling constant J
    is 1.
    
    T : float
    Temperature of our material
    
    Used everywhere.
"""
def analytics(T):
    anaE = -8*(np.sinh(8/T))/(3+np.cosh(8/T))
    anaE2 = 64*(np.cosh(8/T)/(np.cosh(8/T)+3))
    anaM = 0
    anaAbsM = (2*np.exp(8/T)+4)/(np.cosh(8/T)+3)
    anaAbsM2 = (8*np.exp(8/T)+8)/(np.cosh(8/T)+3)
    
    Evar = anaE2 - anaE**2
    Mvar =  anaAbsM2 - anaAbsM**2
    
    cV = CV(Evar, T)
    X = chi(Mvar, T)
    
    store = np.array((anaE, anaAbsM, cV, X))
    return store
"""
    Our probability calculation. Returns a histogram plot.
    
    arr : array-like
    Our input array of all different energies
    L : integer
    The dimension of our lattice, L x L
    T : float
    The temperature of our material
    std : integer
    An indicator for where our function stabilized.
    
    Used for part d)
"""
def probability(arr, L, T, std):
    arr /= L**2
    weights=np.ones_like(arr[std:])/len(arr[std:])
    print('Computed Variance:', np.var(arr))
    plt.hist(arr[std:], bins=20, weights=weights, histtype='bar', ec='black', color='r')
    plt.title('Probability of finding E for T=%g\nafter stabilizing at %i' %(T, std))
    plt.xlabel('Energy')
    plt.ylabel('Probability')
    plt.show()
    
"""
    Unused function for caluclating the speedup given multithreading
    
    t1 : float
    Time taken to calculate for one thread
    tn : float
    Time taken to calculate for n threads
    
    Used nowhere
"""
def speedup(t1, tn):
    return t1/tn
"""
    The heat capacity / specific heat function
    
    Evar : float
    The variance of our energy E
    T : float
    The temperature of our material
    
    Called (and returned) when the boolean extras=True in sorter()
    Additionally used for various plotting purposes
"""
def CV(Evar, T):
    return Evar/(T**2)
"""
    The magnetic susceptibility
    
    Mvar : float
    The variance of our magnetization M
    T : float
    The temperature of our material
    
    Called (and returned) when the boolean extras=True in sorter()
    Additionally used for various plotting purposes
"""
def chi(Mvar, T):
    return Mvar/(T)

