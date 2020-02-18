#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 10:26:14 2020

@author: Sean Dillon

Symplectic Integrator Test

This will (hopefully) successfully test our new symplectic integrator code.
    Pseudocode:
    imports
    functions:
        right hand side function:
            6 vector equations expanded into 18 when I realized solve_ivp 
            does not like vector equations
            6 state equations:
                position_1 (array of size 3)
                position_2
                position_3 
                velocity_1
                velocity_2
                velocity_3
            6 first order differential equations representing two coupled 
            second order equations:
                rhs_1 = velocity_1
                rhs_2 = velocity_2
                rhs_3 = velocity_3
                rhs_4 = G*M2*(x2 - x1)/|x2 - x1|^3 + G*M3*(x3 - x1)/|x3 - x1|^3
                rhs_5 = G*M3*(x3 - x2)/|x3 - x2|^3 + G*M1*(x1 - x2)/|x1 - x2|^3
                rhs_6 = G*M2*(x2 - x3)/|x2 - x3|^3 + G*M1*(x1 - x3)/|x1 - x3|^3
    initial parameters/time setup/initial conditions
    Calculate numerical solutions
    Comparison of numerical methods

"""
import numpy as np
import collections
from symplectic_integrator import  __make_ruth4_update_step_coefficients, integrate
from three_body_problem import rhs
import matplotlib.pyplot as plt

G = 4*np.pi**2 # AU^2/(year^2 * Msun)
M = 1 # Solar masses

t_min = 0 # years
t_max = 1 # years
N = 10000 # fixed 100 steps so animation isn't too slow
times = np.linspace(t_min, t_max, N)
init_cond = np.zeros(18)
init_cond[0] = 1 # r1x; Astronomical units
init_cond[1] = -0.5 # r2x
init_cond[2] = -0.5 # r3x
init_cond[3] = 0 + 1 # v1x
init_cond[4] = -3.3*np.sqrt(3) + 1 # v2x
init_cond[5] = 3.3*np.sqrt(3) + 1 # v3x
init_cond[7] = np.sqrt(3)/2 # r2y
init_cond[8] = -np.sqrt(3)/2 # r3y
init_cond[9] = 6.6 # v1y
init_cond[10] = -3.3 # v2y
init_cond[11] = -3.3 # v3y
init_cond[15] = 1 # v1z (setting all to 1 changes the results interestingly)
init_cond[16] = 1 # v2z
init_cond[17] = 1 # v3z


