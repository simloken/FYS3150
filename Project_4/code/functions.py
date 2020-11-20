import numpy as np

def period(i, mx, add):
    return (i+mx+add)%mx

def monteCarlo(T, L, cycles):
    print('Carlo called\nTemp:', T)
    spintrix = np.zeros( (L,L), np.int8) + 1
    
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
                 
    for i in range(cycles):
        print(i)
        for s in range(L**2):
            x = int(np.random.random()*L)
            y = int(np.random.random()*L)
            deltaE = 2*spintrix.item(x,y)*\
                     (spintrix.item(period(x,L,-1), y) +\
                      spintrix.item(period(x,L,1),  y) +\
                      spintrix.item(x, period(y,L,-1)) +\
                      spintrix.item(x, period(y,L,1)))
            if np.random.random() <= w[deltaE+8]:
                spintrix[x,y] *= -1
                M += 2*spintrix[x,y]
                E += deltaE
        #Update expectation values
        Eex    += E
        E2ex   += E**2
        Mex    += M
        M2ex   += M**2
        absMex += int(np.fabs(M))


    #Normalize average values
    Eex       /= float(cycles);
    E2ex      /= float(cycles);
    Mex       /= float(cycles);
    M2ex      /= float(cycles);
    absMex    /= float(cycles);
    #Calculate variance and normalize to per-point and T
    L2 = L*L
    E_variance  = (E2ex-Eex*Eex)/float(L2*T*T);
    M_variance  = (M2ex-Mex*Mex)/float(L2*T);
    #Normalize returned averages to per-point
    Eex       /= float(L2);
    Mex       /= float(L2);
    absMex    /= float(L2);
    
    return (T, Eex, E_variance, Mex, M_variance, absMex)

def TtoCarlo(T1, T2, TStep, L, cycles):
    T = np.arange(T1, T2+TStep, TStep)
    TL = len(T)
    j = 0
    store = np.zeros((TL, 6))
    for i in T:
        print(i)
        store[j,:] = monteCarlo(i, L, cycles)
        j += 1
    return store

    
def normer(Eex, E2ex, Mex, M2ex, absMex, L, cycles):
    Eex       /= float(cycles)
    E2ex      /= float(cycles)
    Mex       /= float(cycles)
    M2ex      /= float(cycles)
    absMex    /= float(cycles)
    L2 = L*L
    E_variance  = (E2ex-Eex*Eex)/float(L2*T*T)
    M_variance  = (M2ex-Mex*Mex)/float(L2*T)
    Eex       /= float(L2)
    Mex       /= float(L2)
    absMex    /= float(L2)
    
    return (Eex, E_variance, Mex, M_variance, absMex)