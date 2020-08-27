#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 14:58:06 2020

@author: leonard
"""

hbar = 1.0545718*10**(-34)           #hbar in Js
B0 = 2.16                             #Magnetic field in polariser in T
gn = 1.83247171*10**8                #absolute value of neutron gyromagnetic ratio in Hz/T
mn = 1.67492749804*10**(-27)         #neutron mass in kg
e = 1.602176634*10**(-19)            #elementary charge in C
neV = e*10**(-9)
VF = 210*neV                         #Iron Fermi Potential in J
T = 193*neV                          #Temperature of the UCNs in Joule
nFe = 7.874*10**3/(56*mn)            #number density of Fe-56 at room temperature in m^-(3) [Davis, J.R. Metals Handbook Desk Edition. ASM, 1998.]
sigma_abs = 259*10**(-30)            #UCN absorbtion X-section of Fe-56 at room temperature in m^2
sigma_tot = 1242*10**(-30)           #UCN total reaction cross section with Fe-56 at room temperature in m^2
