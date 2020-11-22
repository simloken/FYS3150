from functions import monteCarlo, sorter, Stabilize, normer
import multiprocessing as mp
import numpy as np
import itertools as it

import matplotlib.pyplot as plt

if __name__ == '__main__':
    threads = int(input('Please enter the number of threads you would like to use:\n'))
    if threads > mp.cpu_count():
        raise ValueError('Number of Threads cannot exceed available threads.\nMax threads: %i' %(mp.cpu_count))
    T = 1; L = 2; state=False; accept = True
    cycles = [1000, 10000, 100000, 1000000, 10000000]
    noPar = [10, 100]
    nstore = np.zeros((len(noPar), 6))
    store = np.zeros((len(noPar), 6))
    accepted = np.zeros((len(cycles)+len(noPar)))
    for i in range(len(noPar)):
        if accept == False:
            nstore[i] = monteCarlo(T,L,noPar[i], False, state, accept)
        elif accept == True:
            nstore[i], accepted[i] = monteCarlo(T,L,noPar[i], False, state, accept)
        store[i] = sorter(nstore, T, L, noPar[i], False, True)
        
    preRes = np.zeros((len(cycles), threads, 6))
    with mp.Pool(processes=threads) as pool:
        for i in range(len(cycles)):
            nCycles = cycles[i] // threads
            results = [pool.apply_async(monteCarlo, (T,L, nCycles, True, state, accept)) for iter in range(threads)]
            if accept == False:
                preRes[i] = np.array(([res.get() for res in results])) #extracts all values
            elif accept == True:
                a = np.array(([res.get() for res in results])) #extracts all values
                b = np.array(a[:][:,0]) 
                j = 0
                for k in b:
                    preRes[i][j] = k
                    j += 1
                accepted[i+len(noPar)] = np.sum(a[:][:,1])
        newRes = np.zeros((len(cycles)+len(noPar), 6))
        for i in range(len(store)):
            newRes[i] = store[i]
        j = 0
        for i in preRes:
            newRes[j+len(noPar)] = sorter(i, T, L, cycles[j], False, True)
            j += 1
        for i in cycles:
            noPar.append(i)
        cycles = noPar
        Stabilize(newRes, cycles, L, T) #is as likely
        plt.plot(cycles, accepted)
        plt.show()