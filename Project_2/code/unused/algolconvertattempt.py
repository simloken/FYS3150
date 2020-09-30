# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 14:11:52 2020

@author: Simen
"""

"""
def bisection(a,d,n,eps1,rel=0.1*2**(-217)):
    beta = np.zeros(n)
    beta = a*a
    beta[1] = a[1] = 0
    xmin = d[n-1] - abs(a[n-1])
    xmax = d[n-1] + abs(a[n-1])
    for i in range(n-1,1,-1):
        h = abs(a[i-1])+abs(a[i])
        if d[i] + h > xmax:
            xmax = d[i]+h
        elif d[i] - h < xmin:
            xmin = d[i]-h
    if xmin + xmax > 0:
        eps2 = rel*xmax
    else:
        eps2 = rel*-xmin
    if eps1 <= 0:
        eps1 = eps2
    eps2 = 0.5*eps1+7*eps2
    x0 = xmax
    for i in range(m1,m2):
        x[i] = xmax
        wu[i] = xmin
    z = 0
    for k in range(m2, m1, -1):
        xu = xmin
        for i in range(k, m1, -1):
            if xu<wu[i]:
                if xu == wu[i]:
                    break
        if  x0 > x[k]:
            x0 = x[k]
        x1 = range((xu+x0)/2) 
        while x0 - xu>2*rel*(abs(xu)+abs(x0))+ eps1:
            z = z+1
            a = 0; q= 1
            for i in range(1,n):
                if q != 0:
                    bq = beta[i]/q
                else:
                    bq = abs(b[i]/rel)
                q = c[i] - x1- bq
            if a<k:
                if a<m1:
                    xu = wu[m1] = x1
                else:
                    xu = wu[a+1]=x1
                    if x[a]>x1:
                        x[1]=x1
                    
            else:
                x0 = x1
            x[k] = (x0+xu)/2
    return x
"""