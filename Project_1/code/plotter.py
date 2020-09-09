# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
file1 = open('gsData.txt', 'r'); file2 = open('ssData.txt', 'r'); file3 = open('LUData.txt', 'r')
f1 = file1.readlines(); f2 = file2.readlines(); f3 = file3.readlines()
sols = [f1, f2]
data = []; data2 = []; data3 = []
l = 0
for i in sols: 
    for j in i:
        for k in range(1,4):
            if l == 0:
                data.append(j.split(',')[k])
            elif l == 1:
                    data2.append(j.split(',')[k])
    l += 1
ndata = []; ndata2 = []; ndata3 = []
for i in data:
    ndata.append(i.strip())
for i in data2:
    ndata2.append(i.strip())
data = np.zeros((int(len(ndata)/3),3))
data2 = np.zeros((int(len(ndata2)/3),3))
"""
I am so very sorry for this block of code, this program was a last minute
addition so I had to make something really quick, hence why it's so sloppy.
"""
ndatat = ndata[0::3]
k = 0
for i in ndatat:
    data[k][0] = i
    k += 1
ndatat = ndata[1::3]
k = 0
for i in ndatat:
    data[k][1] = i
    k += 1
ndatat = ndata[2::3]
k = 0
for i in ndatat:
    data[k][2] = i
    k += 1
ndatat = ndata2[0::3]
k = 0
for i in ndatat:
    data2[k][0] = i
    k += 1
ndatat = ndata2[1::3]
k = 0
for i in ndatat:
    data2[k][1] = i
    k += 1
ndatat = ndata2[2::3]
k = 0
for i in ndatat:
    data2[k][2] = i
    k += 1
for i in f3:
    for k in range(2):
        data3.append(i.split(",")[k])
for i in data3:
    ndata3.append(i.strip())
data3 = ndata3[1::2]
ndata3 = data3
data3 = np.zeros(len(data3))
k = 0
for i in ndata3:
    data3[k] = i
    k += 1
"""
Plotting
"""
n1 = [1,2,3,4,5,6,7]
plt.plot(n1, data[:,2]); plt.plot(n1, data2[:,2]); plt.plot(n1[:4], data3)
plt.xlabel("n power"); plt.ylabel("Seconds [s]")
plt.legend(["General Solution", "Special Solution", "LU-Decomposition"])
plt.title("Time comparison between different methods for solving linear equations")
plt.show()

plt.figure()

plt.plot(data[:,1], data[:,0])
plt.plot(data2[:,1], data2[:,0])
plt.legend(["General Solution", "Special Solution"])
plt.title('Relative Error given stepsize h')
plt.xlabel('$Log_{10}(h)$'); plt.ylabel('$Log_{10}(Epsilon)$')
plt.show()