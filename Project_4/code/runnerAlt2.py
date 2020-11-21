from functions import TtoCarlo, monteCarlo
import multiprocessing as mp
import numpy as np
from itertools import product
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
    print('\n\n')
    
    T = np.arange(T1, round(T2+TStep), TStep)
    T = np.around(T, 1)
    if threads > len(T): #excessive threads for area of interest
        threads = len(T)
    with mp.Pool(processes=threads) as pool:
        t0 = time.perf_counter()
        results = [pool.apply_async(monteCarlo, (T1,L,cycles)) for i in range(threads)]
        print([res.get() for res in results])
        print('\nTime:\n',time.perf_counter() - t0)
        print('\nComplete!')
        last = (input('Press Enter to close program'))