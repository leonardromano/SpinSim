#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 12:43:29 2020

@author: leonard
"""
import numpy as np

def zero(x):
    return 0

def constant(x):
    return 1

def linear(x):
    return 0.1 + 0.9*x

def sigmoid(x):
    return 1/(1 + (10-1)*np.exp(-10*x))

def exp(x):
    return 10**(3*(x-1))