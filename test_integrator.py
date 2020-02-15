#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 10:26:14 2020

@author: Sean Dillon

Symplectic Integrator Test

This will (hopefully) successfully test our new symplectic integrator code.
"""

import numpy as np
import collections
from symplectic_integrator import  __make_ruth4_update_step_coefficients, integrate


G = 4*np.pi**2 # AU^2/(year^2 * Msun)
M = 1 # Solar masses

t_min = 0 # years
t_max = 1 # years
N = 10000 # fixed 100 steps so animation isn't too slow
times = np.linspace(t_min, t_max, N)
