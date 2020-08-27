#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 12:09:21 2020

@author: leonard
"""

import numpy as np
import Integration_Utility as IU

def Euler(N, u, L, B, w, S0, out = False):
    "Integrates the Propagator equation dP=MPds and returns the result"
    Prop = np.identity(3)
    if np.all(S0==0):
        return Prop
    else:
        if out == False:
            for i in range(N):
                Prop += np.dot(IU.M(B, i/N, w, L, u, np.dot(Prop, S0)[2], N),Prop)*(1/N)
        else:
            for i in range(N):
                Prop += np.dot(IU.M(B, 1-i/N, w, L, u, np.dot(Prop, S0)[2], N),Prop)*(1/N)
        return Prop

def RK4(N, u, L, B, w, S0, out = False):
    "Integrates the Propagator equation using RK4 dP=MPds and returns the result"
    Prop = np.identity(3)
    h = 1/N
    if np.all(S0==0):
        return Prop
    elif L==0:
        return np.diag([-1,1,-1])
    else:
        if out == False:
            for i in range(N):
                k1 = np.dot(IU.M(B, i*h, w, L, u, np.dot(Prop, S0)[2], N),Prop)
                k2 = np.dot(IU.M(B, (i+1/2)*h, w, L, u, np.dot(Prop+h*k1/2, S0)[2], N),Prop+k1*h/2)
                k3 = np.dot(IU.M(B, (i+1/2)*h, w, L, u, np.dot(Prop+h*k2/2, S0)[2], N),Prop+k2*h/2)
                k4 = np.dot(IU.M(B, (i+1)*h, w, L, u, np.dot(Prop+h*k3, S0)[2], N),Prop+h*k3)
                Prop += h*(k1 + 2*k2+2*k3+k4)/6
        else:
           for i in range(N):
                k1 = np.dot(IU.M(B, 1-i*h, w, L, u, np.dot(Prop, S0)[2], N),Prop)
                k2 = np.dot(IU.M(B, 1-(i+1/2)*h, w, L, u, np.dot(Prop+h*k1/2, S0)[2], N),Prop+k1*h/2)
                k3 = np.dot(IU.M(B, 1-(i+1/2)*h, w, L, u, np.dot(Prop+h*k2/2, S0)[2], N),Prop+k2*h/2)
                k4 = np.dot(IU.M(B, 1-(i+1)*h, w, L, u, np.dot(Prop+h*k3, S0)[2], N),Prop+h*k3)
                Prop += h*(k1 + 2*k2+2*k3+k4)/6
        return Prop

def RK4EulerAngles(N, u, L, B, w, S0, out = False):
    alpha = np.array([0.,0.,0.])
    h = 1/N
    if np.all(S0==0):
        alpha*=0.
    elif L==0:
        alpha += np.array([0.,np.pi,0.])
    else:
        if out == False:
            alpha += h*np.array([IU.kappa(0, L, u, B, S0[2], N) - IU.dw(w), \
                                 IU.dw(w), IU.dw(w)])
            for i in range(1,N):
                k1 = IU.dAlpha(i*h, alpha, L, u, B, S0, w, N)
                k2 = IU.dAlpha((i+1/2)*h, alpha + h*k1/2, L, u, B, S0, w, N)
                k3 = IU.dAlpha((i+1/2)*h, alpha + h*k2/2, L, u, B, S0, w, N)
                k4 = IU.dAlpha((i+1)*h, alpha + h*k3, L, u, B, S0, w, N)
                alpha += h*(k1 + 2*k2+2*k3+k4)/6
        else:
            alpha += h*np.array([IU.kappa(1, L, u, B, S0[2], N) - IU.dw(w), \
                                 IU.dw(w), IU.dw(w)])
            for i in range(1,N):
                k1 = IU.dAlpha(1-i*h, alpha, L, u, B, S0, w, N)
                k2 = IU.dAlpha(1-(i+1/2)*h, alpha + h*k1/2, L, u, B, S0, w, N)
                k3 = IU.dAlpha(1-(i+1/2)*h, alpha + h*k2/2, L, u, B, S0, w, N)
                k4 = IU.dAlpha(1-(i+1)*h, alpha + h*k3, L, u, B, S0, w, N)
                alpha += h*(k1 + 2*k2+2*k3+k4)/6
    return IU.Phi(*alpha)

def chooseIntegrator(Integ):
    if Integ == "Euler":
        return Euler
    elif Integ == "RK4":
        return RK4
    elif Integ == "AngleRK4":
        return RK4EulerAngles
    else:
        print("Unknown argument: Integrator")