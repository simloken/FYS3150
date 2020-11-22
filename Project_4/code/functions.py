import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
def period(i, mx, add):
    return (i+mx+add)%mx

def monteCarlo(T, L, cycles, progress, state=False, accept=False):
    if progress==True:
        print('ID %i working through %i cycles' %(int(mp.current_process()._identity[0]), cycles))
    if state == False:
        spintrix = np.random.choice([-1,1], size=(L,L)) #either up or down
    else:
        spintrix = np.full((L,L), fill_value=1)
        
    nonTen = False
    count = 0
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
    if progress==True:
        print('===|ID %i finished|===' %(int(mp.current_process()._identity[0])))
    if accept == True:
        return store, count
    else:
        return store

def sorter(arr,T, L, cycles, extras, alt=False):
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
                X[i] = chi(sarr[i,4], T[i])
            
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
    M_variance  = (M2ex-Mex**2)/float(L2*T)
    
    Eex       /= float(L2)
    Mex       /= float(L2)
    absMex    /= float(L2)
    
    narr = np.array((T, Eex, E_variance, Mex, M_variance, absMex))
    return narr


def Stabilize(arr, cL, L, T):
    analytical = analytics(T)/(L**2)
    print(arr)
    plt.figure()
    plt.semilogx(cL, analytical[0]*np.ones(len(cL)),'r.-', label='Analytic Mean E')
    plt.semilogx(cL, arr[:,1], 'gx-', label='Num. Mean E')
    plt.legend()
    plt.show()
    plt.figure()
    plt.semilogx(cL, analytical[1]*np.ones(len(cL)), 'r.-', label='Analytic Mean M')
    plt.semilogx(cL, arr[:,5], 'gx-', label='Num.Mean M')
    plt.legend()
    plt.show()
    plt.figure()
    plt.semilogx(cL, analytical[2]*np.ones(len(cL)), 'r.-', label='Analytic Specific Heat')
    plt.semilogx(cL, CV(arr[:,2], T), 'gx-', label='Num. Specific Heat')
    plt.legend()
    plt.show()
    plt.figure()
    plt.semilogx(cL, analytical[3]*np.ones(len(cL)),'r.-', label='Analytic Mag. Susc.')
    plt.semilogx(cL, chi(arr[:,4], T), 'gx-', label='Num. Mag. Susceptibility')
    plt.legend()
    plt.show()


def analytics(T):
    anaE = -8*(np.sinh(8/T))/(3+np.cosh(8/T))
    anaE2 = -64*(np.cosh(8/T)/(np.cosh(8/T)+3))
    anaM = 0
    anaAbsM = (2*np.exp(8/T)+4)/(np.cosh(8/T)+3)
    anaAbsM2 = (8*np.exp(8/T)+8)/(np.cosh(8/T)+3)
    Evar = anaE2 - anaE**2
    Mvar =  anaAbsM**2 - anaAbsM2
    
    cV = CV(Evar, T)
    X = chi(Mvar, T)
    store = np.array((anaE, anaAbsM, cV, X))
    return store

def speedup(t1, tn):
    return t1/tn

def CV(Evar, T):
    print(Evar)
    return Evar/(T**2)

def chi(Mvar, T):
    return Mvar/(T)
