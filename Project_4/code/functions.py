import numpy as np
import multiprocessing as mp
def period(i, mx, add):
    return (i+mx+add)%mx

def monteCarlo(T, L, cycles, progress):
    print('ID %i working through %i cycles' %(int(mp.current_process()._identity[0]), cycles))
    spintrix = np.zeros( (L,L), np.int8) + 1
    nonTen = False
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
        
        if progress==True:
            if nonTen == True:
                if i % int(cycles/10) == 0 and i != 0: #special case so we don't get percentages higher than 100
                    k += 1
                    print('ID %i at %i%%' %(int(mp.current_process()._identity[0]), int(k*10)))
            else:
                if i % int(cycles/10) == 0:
                    k += 1
                    print('ID %i at %i%%' %(int(mp.current_process()._identity[0]), int(k*10)))

    store = np.array((T, Eex, E2ex, Mex, M2ex, absMex))
    print('===|ID %i finished|===' %(int(mp.current_process()._identity[0])))
    return store

def sorter(arr,T, L, cycles):
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
    return sarr
    
def normer(arr, L, cycles):
    T = arr[0]
    Eex = arr[1]; E2ex = arr[2]; Mex = arr[3]; M2ex = arr[4]; absMex = arr[5]
    
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
    
    narr = np.array((T, Eex, E_variance, Mex, M_variance, absMex))
    
    return narr


def TtoCarlo(T1, T2, TStep, L, cycles):
    T = np.arange(T1, T2+TStep, TStep)
    TL = len(T)
    j = 0
    store = np.zeros((TL, 6))
    for i in T:
        store[j,:] = monteCarlo(i, L, cycles)
        j += 1
    return store

