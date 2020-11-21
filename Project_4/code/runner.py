from functions import TtoCarlo, monteCarlo, sorter
import multiprocessing as mp
import numpy as np
import itertools as it
import time
import sys
import os

    
    
    
if __name__ == '__main__':
    if len(sys.argv) != 7:
        print('Some initial inputs missing!\nRequired Inputs:[threads, L, cycles, T1, T2, TStep]\nPlease re-enter inputs below as requested')
        threads = int(input('Please enter the number of threads you would like to use:\n'))
        if threads > mp.cpu_count():
            raise ValueError('Number of Threads cannot exceed available threads.\nMax threads: %i' %(mp.cpu_count))
        L = int(input('Please enter an integer value for L:\n'))
        cycles = int(input('Please enter an integer value for number of Montecarlo Cycles:\n'))
        T1 = float(input('Intial Temperature:\n'))
        T2 = float(input('Final Temperature:\n'))
        TStep = float(input('Temperature Step Size:\n'))
    else:
        threads = int(sys.argv[1])
        if threads > mp.cpu_count():
            raise ValueError('Number of Threads cannot exceed available threads.\nMax threads: %i' %(mp.cpu_count))
        L = int(sys.argv[2])
        cycles = int(sys.argv[3])
        T1 = float(sys.argv[4])
        T2 = float(sys.argv[5])
        TStep = float(sys.argv[6])
    excess = False; special = False; leftover = False
    dd = 0
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
        if excess != True: #if there is not an excess in threads in relation to temperatures, assign a thread to each temperature
            results = [pool.apply_async(monteCarlo, (temp,L,cycles)) for temp in T]
        elif excess == True: #if there is excess, assign additional threads to chunks/batches of cycles to speed up process
            if special == True:
                    results = ([pool.apply_async(monteCarlo, (temp,L,nCycles)) for temp in T for iter in range(d)])
            else:
                    leftover = True
                    results = ([pool.apply_async(monteCarlo, (temp,L,nCycles)) for temp in T[:len(T)-1] for iter in range(d)])
                    leftover_res = ([pool.apply_async(monteCarlo, (T[len(T)-1],L,nnCycles)) for iter in range(dd)])
                    for i in leftover_res:
                        results.append(i)
                
        preRes = np.array(([res.get() for res in results])) #extracts all values
        t1 = time.perf_counter() - t0
        sortRes = sorter(preRes, T, L, cycles) #sorts and scales/normalizes the values
        print(sortRes)
        print('\nTime:\n%f' %(t1))
        print('\nComplete!')
        cwd = os.getcwd()
        os.chdir("../data/")
        if os.path.isfile('data.txt') == False: #creates files and labels for first run
            file = open('data.txt', 'w+')
            file.write('Cycles\tNSpins\tTemp\tMeanE\t\tEnergyVar\tMeanM\t\tMVar\t\tMeanAbsM\tTime\n')
            file.close()
        file = open('data.txt', 'a')
        for i in sortRes:
            file.write("%i\t%i\t%0.3f\t%f\t%f\t%f\t%f\t%f\t%g\n" %(cycles,L,i[0],i[1],i[2],i[3],
                                                                i[4], i[5], t1))
        file.close()
        os.chdir(cwd)
        last = (input('Press Enter to close program'))