import matplotlib.pyplot as plt
from scipy import optimize
import numpy as np
"""
    My own file-reading function, hence all the extra unecessary (in this case) features
"""
def reader(filename, p=None, d=None):
    if isinstance(filename, str) == False:
        raise ValueError('Filename must be string')
    file = open(filename, 'r')
    f = file.readlines()
    file.close()
    
    if p == None:    
        p = input(str('Does the file use tabs or newlines? \nIf so, input Y. \nIf not, please input N for commas\n'))
        p = p.lower()
        if p not in ['y','n']:
            raise ValueError('Input must be either Y, N or empty!')
            
    if d == None:
        d = int(input('How many variables are in the data set?\n'))
    
    data = []
    if p == 'y':
        for i in f:
                data.append(i.split('\t'))
        f = []  
        getlen = len(data)
        for i in range(len(data)):
            for j in data[i]:
                f.append(j.strip())
        data = np.zeros((d,getlen))
        for i in range(d):
            data[i] = f[i::d]
            
    elif p == 'n':
        for i in f:
            data.append(i.split(','))
        f = []
        for i in data[0]:
            f.append(i.strip())
        data = [float(i) for i in f]
        getlen= len(data)
        f = data
        data = np.zeros((d,int(getlen/d)))
        for i in range(d):
            data[i] = f[i::d]
    
      
    
    data = np.transpose(data) #looks nicer
    return data

data = reader('../data/timeData.txt', 'y', 9)
threads = [1,2,4,8,16]
times = data[:,8]
times = times[0::2]
def func(x,a,b, c): #looks a bit like this function
    return a*1/x*b +c
params, _ = optimize.curve_fit(func, threads, times)
plt.scatter(threads, times, color='red', label='Data point')
plt.plot(threads, func(threads, params[0], params[1], params[2]), label='Fit line')
plt.xlabel('CPU Threads')
plt.ylabel('Time [s]')
plt.title('Time spent calculating a for L=20, Cycles = 100000\nfor differing number of active CPU Threads')
plt.legend()
plt.show()