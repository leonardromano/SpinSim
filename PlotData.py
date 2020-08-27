#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 14:38:45 2020

@author: leonard
"""

import matplotlib.pyplot as plt
import numpy as np

path = "<path-to-data>"
files = ["flipper on RK4" , "flipper off RK4"]
title = "L=50nm, B=const., w=5"

for file in files:
    x, y, N, t = np.loadtxt(path + file, unpack=True)
    plt.plot(x,y, label=file[:-4])
plt.xlabel("speed [m/s]")
plt.ylabel("Polarisation")
plt.title(title)
plt.legend()
plt.xlim(3.9,12)
plt.ylim(-1.05, 1.05)
plt.savefig(path+"P-u.png")
plt.show()

for file in files:
    x, y, N, t = np.loadtxt(path + file, unpack=True)
    plt.plot(x,N, label=file[:-4])
plt.xlabel("speed [m/s]")
plt.ylabel("dN/du")
plt.title(title)
plt.legend()
plt.xlim(3.9,12)
#plt.ylim(-1.05, 1.05)
plt.savefig(path+"N-u.png")
plt.show()

for file in files:
    x, y, N, t = np.loadtxt(path + file, unpack=True)
    plt.plot(x, t, label=file[:-4])
plt.xlabel("speed [m/s]")
plt.ylabel("N/N0")
plt.title(title)
plt.legend()
plt.xlim(3.9,12)
plt.ylim(0, 1.01)
plt.savefig(path+"T-u.png")
plt.show()


