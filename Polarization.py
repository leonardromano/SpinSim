#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 14:13:10 2020

@author: leonard
"""
import numpy as np
from numpy import sin
from numpy import cos
import matplotlib.pyplot as plt

def dk(k):
    return (2*k+1)*np.pi
    
def z(x,k):
    return np.sqrt(dk(k)**2 + x**2)

def sinc(x):
    return np.sinc(x/np.pi)

def P(x, l, k):
    zeta = z(x,k)
    return (dk(k)**2)*x*sin(l*x)*sinc(zeta)*sinc(zeta/2)**2 + \
        (dk(k)**2)*cos(l*x)*((x/2)**2*sinc(zeta/2)**4 - sinc(zeta)**2) + \
            ((dk(k)**2*cos(zeta) + x**2)/zeta**2)**2

l = 5
x = np.arange(0, 100, 0.01)
ks = [0,5]
for k in ks:
    plt.plot(x, P(x, l, k), label = "k = %i"%(k))
plt.xlabel(r"$\Delta\varphi$")
plt.ylabel("P")
plt.legend(loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=3)
plt.ylim(-1,1)    
plt.show()

ls = [2, 100]
for l in ls:
    plt.plot(x, P(x, l, 10), label = "l = %2.1f"%(l))
    plt.xlabel(r"$\Delta\varphi$")
    plt.ylabel("P")
    plt.legend(loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=3)
    plt.ylim(-1,1)    
    plt.show()
        