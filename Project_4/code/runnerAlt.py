from functions import TtoCarlo
import multiprocessing as mp
from itertools import product
import time
import sys

#divides total cycles for amount of threads. Not efficient and probably not correct
#comment further later
    
    
    
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
    nCycles = cycles // threads
    nCycles = int(nCycles)
    t0 = time.perf_counter()
    lst = []
    
    with mp.Pool(processes=threads) as pool:
        results = [pool.apply_async(TtoCarlo, (T1, T2, TStep, L, cycles,)) for i in range(threads)]
        print([res.get() for res in results])
        print('\nTime:\n',time.perf_counter() - t0)
        print('\nComplete!')
        last = (input('Press Enter to close program'))