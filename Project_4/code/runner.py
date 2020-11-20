from functions import TtoCarlo, monteCarlo
import multiprocessing as mp
import numpy as np
import itertools as it
import time

    
    
    
if __name__ == '__main__':
    threads = int(input('Please enter the number of threads you would like to use:\n'))
    if threads > mp.cpu_count():
        raise ValueError('Number of Threads cannot exceed available threads.\nMax threads: %i' %(mp.cpu_count))
    
    L = int(input('Please enter an integer value for L:\n'))
    cycles = int(input('Please enter an integer value for number of Montecarlo Cycles:\n'))
    T1 = float(input('Intial Temperature:\n'))
    T2 = float(input('End Temperature:\n'))
    TStep = float(input('Temperature Step:\n'))
    excess = False; special = False
    
    T = np.arange(T1, T2+TStep, TStep)
    np.round(T,1)
    if threads > len(T): #if threads are greater than temperature then
        excess = True #use excess threads to split up cycles in smaller chunks and mend later
        if threads % len(T) == 0: #special case where we get a neater solution
            special = True
            d = threads//len(T)
            nCycles = cycles // d
        else:
            d = int(np.floor(threads/len(T)))
            nCycles = cycles // d
            dd = int(threads - len(T)*d)
            nnCycles = cycles // dd
            
    with mp.Pool(processes=threads) as pool:
        t0 = time.perf_counter()
        if excess != True:
            results = [pool.apply_async(monteCarlo, (temp,L,cycles)) for temp in T]
        elif excess == True:
            if special == True:
                    results = ([pool.apply_async(monteCarlo, (temp,L,nCycles)) for temp in T for iter in range(d)])
            else:
                    results = ([pool.apply_async(monteCarlo, (temp,L,nCycles)) for temp in T[:len(T)-1] for iter in range(d)])
                    leftover_res = ([pool.apply_async(monteCarlo, (T[len(T)-1],L,nnCycles)) for iter in range(dd)])
                    for i in leftover_res:
                        results.append(i)
                    
                
        print([res.get() for res in results])
        print('\nTime:\n%f' %(time.perf_counter() - t0))
        print('\nComplete!')
        last = (input('Press Enter to close program'))