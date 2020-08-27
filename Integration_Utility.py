#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 13:19:11 2020

@author: leonard
"""
import numpy as np
from numpy import sin,cos,tan,cosh,sinh, exp
from Constants import *

#useful functions for the quantum transmission function
def transmission(u,d, s):
    "calculates the transmission parameter as determined from a 1D-potential wall"
    E = mn*u**2/2
    if s == "up":
        V = VF+hbar*gn*B0/2
    elif s=="down":
        V = VF-hbar*gn*B0/2
    else:
        print("input not defined! s can only be 'up' or 'down'!")
        V=0
    x = V/E
    if x==1:
        return 1
    elif V > E:
        k = np.sqrt(2*mn*(V-E))/hbar
        return (x-1)/((x-1)*cosh(k*d)**2 + (1-x/2)**2 * sinh(k*d)**2)
    else:
        k=np.sqrt(2*mn*(E-V))/hbar
        return (1-x)/((1-x)*cos(k*d)**2 + (1-x/2)**2 * sin(k*d)**2)
    
def polariser(u, d, P, Nup, Ndown):
    "updates the polarisation according to the polariser formula"
    tup = transmission(u,d, "up")
    tdown = transmission(u,d, "down")
    if abs(P)==1:
        R=P
    else:
        R = (tup-tdown + (tup + tdown)*P)/(tup+tdown+(tup-tdown)*P)
    if R>1:
        return 1, 0, Ndown
    elif R<-1:
        return -1, Nup, 0
    else:
        return R, tdown*Nup, tup*Ndown

def polarise(S, u, d, Nup, Ndown):
    "Takes the incoming polarisation vector and applies the polariser to it"
    P0 = -S[2]          # The "-" comes from the fact that the polariser is 
                        #defined wrt to the unflipped coordinate system
    P, Nup, Ndown = polariser(u,d, P0, Nup, Ndown)
    if abs(P0) == 1:
        S[0] = 0
        S[1] = 0
    else:
        S[0] = np.sqrt((1-P**2)/(1-P0**2))*S[0]
        S[1] = np.sqrt((1-P**2)/(1-P0**2))*S[1]
    S[2] = -P
    return S, Nup, Ndown

#Free precession
def rotationz(phi):
    "A rotation around the z-axis with angle phi"
    return np.array([[cos(phi), sin(phi), 0],\
                     [-sin(phi), cos(phi), 0],\
                         [0,0,1]])
        
def constantOne(x):
    return 1

def k0(P, u):
    "The flipping frequency inside the polariser"
    return gn*B0/u/v(1, u, constantOne, P)

#useful functions for the RK4 numerical integration
def v(x, u, B, P):
    "returns the velocity at position x in m/s"
    y = 1 + (hbar*gn*B0*B(x)*P- 2*VF*np.heaviside(x-0.5, 1))/mn/u**2
    if y <=0:
        return np.sqrt(abs(y))
    return np.sqrt(y)

def dw(w):
    "Returns the differential of the winding phase of the magnetic field"
    return (2*w+1)*np.pi

def kappa(x, L, u, B, P, N):
    "Returns the differential of the flipping phase"
    if L==0:
        return 0
    y = (gn*B0*L)*B(x)/v(x,u,B,P)/u
    if y >= N:
        return N
    else:
        return y

def M(B, x, w, L, u, P, N):
    "Returns the differential of the Matrix propagator"
    if P>1:
        P=1
    elif P<-1:
        P=-1
    k = kappa(x, L, u, B, P, N)
    return  np.array([[0,k, -dw(w)],\
                      [-k, 0, 0],\
                          [dw(w), 0, 0]])

#Useful functions for the particle count part
def maxwellBoltzmann(u):
    "Maxwell-Boltzmann Distribution at temperature T and velocity u"
    return 4*np.pi*(mn/2/np.pi/T)**(3/2) * u**2 *exp(-mn*u**2 /2/T)

def updateNs(N, P):
    "# of spins up and down for total number N and polarisation P"
    return N/2 * (1+P), N/2 * (1-P)

## Useful functions for the classical transmission
def transmissionCoefficients(u, d, pol):
    "calculates the transmission and reflection coefficients due to scattering and absorbtion"
    if pol == "up":
        y = 1 + (hbar*gn*B0 - 2*VF)/mn/u**2
    elif pol == "down":
        y = y = 1 - (hbar*gn*B0 + 2*VF)/mn/u**2
    else:
        y=0
    if y <=0:
        return 0, 0, 1
    else:
        v=u*np.sqrt(y)
        return exp(-nFe*sigma_abs*d*2200/v), 4*u*v/(u+v)**2,((u-v)/(u+v))**2

def absorb(S, u, d, Nup, Ndown):
    "Updates the number of neutrons and their polarisation after absorption"
    xup, Tup, Rup = transmissionCoefficients(u, d, "up")
    xdown, Tdown, Rdown = transmissionCoefficients(u, d, "down")
    alpha = np.exp(-nFe*sigma_tot*d)
    Nup*=alpha*xup*Tup**2/(1-(alpha*xup*Rup)**2)
    Ndown*=alpha*xdown*Tdown**2/(1-(alpha*xdown*Rdown)**2)
    P = (Nup-Ndown)/(Nup+Ndown)
    if abs(S[2]) == 1:
        S[0] = 0
        S[1] = 0
    else:
        S[0] = np.sqrt((1-P**2)/(1-S[2]**2))*S[0]
        S[1] = np.sqrt((1-P**2)/(1-S[2]**2))*S[1]
    S[2] = P
    return S, Nup, Ndown

# Useful functions for the numerical integration with Euler angles

def Phi(alpha, beta, gamma):
    return np.array([[cos(alpha+gamma) + (cos(beta) - 1)*cos(alpha)*cos(gamma), \
                      sin(alpha+gamma) + (cos(beta) - 1)*cos(alpha)*sin(gamma), \
                          -sin(beta)*cos(alpha)],\
                     [-sin(alpha+gamma) - (cos(beta) - 1)*sin(alpha)*cos(gamma), \
                      cos(alpha+gamma) - (cos(beta) - 1)*sin(alpha)*sin(gamma), \
                          sin(beta)*sin(alpha)],\
                     [sin(beta)*cos(gamma), sin(beta)*sin(gamma), cos(beta)]]) 
                

def dAlpha(x, alpha, L, u, B, S0, w, N):
    return np.array([kappa(x, L, u, B, np.dot(Phi(*alpha),S0)[2],N) - sin(alpha[0])/tan(alpha[1]) * dw(w),\
                     cos(alpha[0])*dw(w), \
                     sin(alpha[0])/sin(alpha[1]) *dw(w)])