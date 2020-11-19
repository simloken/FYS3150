from functions import TtoCarlo, monteCarlo
import multiprocessing as mp
from itertools import product
import time
import sys
if len(sys.argv) > 1:
    if len(sys.argv) == 2:
        threads = int(sys.argv[1])
        print('Threads: %i' %(sys.argv[1]))
        L = int(input('Please enter an integer value for L:\n'))
        cycles = int(input('Please enter an integer value for number of Montecarlo Cycles:\n'))
    elif len(sys.argv) == 3:
        threads = int(sys.argv[1])
        L = int(sys.argv[2])
        print('Threads: %i' %(sys.argv[1]))
        print('L: %i' %(sys.argv[2]))
        cycles = int(input('Please enter an integer value for number of Montecarlo Cycles:\n'))
    elif len(sys.argv) == 4:
        threads = int(sys.argv[1]); L = int(sys.argv[2]); cycles = int(sys.argv[3])
        print('Threads: %i\nL: %i\nMontecarlo Cycles: %i' %(threads,L,cycles))
    threads = int(sys.argv[1]); L = int(sys.argv[2]); cycles = int(sys.argv[3])
else:
    threads = int(input('Please enter the number of threads you would like to use:\n'))
    L = int(input('Please enter an integer value for L:\n'))
    cycles = int(input('Please enter an integer value for number of Montecarlo Cycles:\n'))
    
    
    
if __name__ == '__main__':
    if threads > mp.cpu_count():
        raise ValueError('Number of Threads cannot exceed available threads.\nMax threads: %i' %(mp.cpu_count))
        
    T1 = float(input('Intial Temperature:\n'))
    T2 = float(input('End Temperature:\n'))
    TStep = float(input('Temperature Step:\n'))
    nCycles = cycles // threads
    nCycles = int(nCycles)
    t0 = time.perf_counter()
    with mp.Pool(processes=threads) as pool:
        results = [pool.apply_async(TtoCarlo, (T1, T2, TStep, L, nCycles, )) for i in range(threads)]
    #pool.close()
    #pool.join()
    print([res.get() for res in results])
    print(time.perf_counter() - t0)
    print('Complete')