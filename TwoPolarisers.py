#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 14:43:15 2020

@author: leonard
"""

import numpy as np
import NumericalIntegration as NI
import Integration_Utility as IU

def onePolariser(S, B, u, L, w, N, d, Nup, Ndown, tMod, Integration):
    "updates S by one polariser"
    #1 entering the medium
    S = np.dot(Integration(N, u, L, B, w, S), S)
    Nup, Ndown = IU.updateNs(Nup+Ndown, S[2])
    #2 absorption or QM-transmission
    if tMod == "classical":
        S, Nup, Ndown = IU.absorb(S, u, d, Nup, Ndown)
    elif tMod == "quantum":
        S, Nup, Ndown = IU.polarise(S, u, d, Nup, Ndown)
    else:
        print("Invalid choice for transmission Model. \
              Choose either 'classical' or 'quantum'")
    #4 free precession
    S = np.dot(IU.rotationz(IU.k0(S[2],u)*d), S)
    
    #5 leaving the medium
    S = np.dot(Integration(N, u, L, B, w, S, True),S)
    Nup, Ndown = IU.updateNs(Nup+Ndown, S[2])
    return S, Nup, Ndown

def twoPolarisers(P0, N0, B, u, L, w, N, d, flip=False, \
                  tMod="classical", Integrator="RK4"):
    "P0 through the whole polariser system and returns the Polarisation"
    Integration = NI.chooseIntegrator(Integrator)
    #initialise the particle number according to Maxwell-Boltzmann distribution
    n0 = N0*IU.maxwellBoltzmann(u)
    #initialise polarisation vector
    S0 = np.array([0,0,P0])
    #initialise number of spins up and down
    Nup, Ndown = IU.updateNs(n0, P0)
    #start propagation through system
    S, Nup, Ndown = onePolariser(S0, B, u, L, w, N, d, Nup, Ndown, \
                                 tMod, Integration)
    if flip ==True:
        S[2] = -S[2]
        Nup, Ndown = IU.updateNs(Nup+Ndown, S[2])
    S, Nup, Ndown = onePolariser(S, B, u, L, w, N, d, Nup, Ndown, \
                                 tMod, Integration)
    return S[2], Nup+Ndown, (Nup+Ndown)/n0


