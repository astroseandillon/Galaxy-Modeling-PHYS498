#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 09:29:22 2020

@author: Sean Dillon

Adding a new integration method using symplectic integration as cited in 
Donnely and Rogers, 2005
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import astropy as ap


#Each point has a set of points (q, p), where q = position and 
#p = momentum. we then create a 2n-dimensional set of points

n = 100

q = 'q'
p = 'p'
arr = np.array((1,n))
for i in range(n):
    arr((0, i)) = q + str(i)
    arr((1, i)) = p + str(i)


