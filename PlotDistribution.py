#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 09:44:13 2020

@author: leonard
"""

import matplotlib.pyplot as plt
import numpy as np

mn = 1.67492749804*10**(-27)         #neutron mass in kg
e = 1.602176634*10**(-19) 
neV = e*10**(-9)
T1 = 80*neV
T = 420*neV   

def maxwellBoltzmann(u):
    "Maxwell-Boltzmann Distribution at temperature T and velocity u"
    return 4*np.pi*(mn/2/np.pi/T)**(3/2) * u**2 *np.exp(-mn*u**2 /2/T)

def gauß(u):
    return 1/np.sqrt(2*np.pi*0.2138) * np.exp(-(u-6.08)**2 /2/0.2138)

N = 50000
u = np.arange(0,12, 0.01)
dNdu = N*gauß(u)

plt.plot(u,dNdu)
plt.xlabel("speed [m/s]")
plt.ylabel("dN/du")
plt.title("Initial distribution, T=193neV")
plt.savefig("dNdu-u.png")
plt.show()
