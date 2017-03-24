#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 10:52:01 2016

@author: philipp
"""
from __future__ import division # floating point division by default
import numpy

def gini(data):
    N = len(data)  
    x = [(i+1)/N for i in range(N)]
    data_sorted = sorted(data)
    S = sum(data_sorted)
    data_normalized = [data_sorted[k]/S for k in range(N)]
    y = list(numpy.cumsum(data_normalized))
    # Calculate Gini coefficient
    d_x = x[0]
    y.insert(0,0)
    x.insert(0,0)
    Y = [(y[k]+y[k+1])/2*d_x for k in range(N)]
    I = sum(Y)
    G = 1 - 2*I
    return G, x, y
