#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 13:27:56 2020

@author: leonard
"""

import functions as func
import numpy as np
from TwoPolarisers import twoPolarisers as run
import time
nm = 10**(-9)

def main(P0, N0, us, L, B, d, w, N, tMod, Integrator):
    "Integrates the Equations of Motion and saves the result in a file"
    starttime = time.time()
    OutputTrue = []
    OutputFalse = []
    for u in us:
        Bools = [False, True]
        for b in Bools:
            if b == True:
                OutputTrue.append([u,*run(P0, N0, B, u, L, w, N, d, b, tMod, Integrator)])
            else:
                OutputFalse.append([u,*run(P0, N0, B, u, L, w, N, d, b, tMod, Integrator)])
    
    np.savetxt("Results/flipper on RK4", np.asarray(OutputTrue), delimiter = " ")
    np.savetxt("Results/flipper off RK4", np.asarray(OutputFalse), delimiter = " ")
    endtime = time.time()
    print("Finished in " + str((endtime-starttime)/60) + " mins\n\n")
    
P0 = 0                              #initial Polarisation
N0 = 50000                          #initial number of UCNs
L = 10*nm                           #Length of flipping area -> implement variable length
Bshape = func.constant              #shape of magnetic field (see functions for options)
d = 300*nm                          #Thickness of polariser
w = 0                               #Winding number
N = 100                             #Number of integration steps
tMod = "classical"                  #Transmission model ("classical"--smooth, or "quantum")
Integrator = "AngleRK4"             #Integration method ("Euler", "RK4", "AngleRK4")
# intialize speed range u in [3.91, 12] (different granularities)
u1 = np.arange(3.91, 4.1, 0.01)
u2 = np.arange(4.1, 7.9, 0.01)
u3 = np.arange(7.9, 8.2, 0.01)
u4 = np.arange(8.2, 12, 0.01)
us = np.append(u1, [*u2, *u3, *u4])


main(P0, N0, us, L, Bshape, d, w, N, tMod, Integrator)


        
