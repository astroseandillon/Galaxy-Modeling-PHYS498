#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 10:26:14 2020

@author: Sean Dillon

Symplectic Integrator Test


"""
import numpy as np
import collections
from symplectic_integrator import  __make_ruth4_update_step_coefficients, integrate
import matplotlib.pyplot as plt

plt.close('all')

def f(x):
    return x**4+3*x**4 + 5

G = 4*np.pi**2 # AU^2/(year^2 * Msun)
M1 = 1 # Solar masses
M2 = 100
M3 = 1

# differential equations set up
def rhs(time, state):
    
    # 18 equations: 3 2nd order * 2 for separating, * 3 for vector components
    position_1x = state[0]
    position_2x = state[1]
    position_3x = state[2]
    velocity_1x = state[3]
    velocity_2x = state[4]
    velocity_3x = state[5]
    position_1y = state[6]
    position_2y = state[7]
    position_3y = state[8]
    velocity_1y = state[9]
    velocity_2y = state[10]
    velocity_3y = state[11]
    position_1z = state[12]
    position_2z = state[13]
    position_3z = state[14]
    velocity_1z = state[15]
    velocity_2z = state[16]
    velocity_3z = state[17]
    
    position_1 = np.array([position_1x, position_1y, position_1z])
    position_2 = np.array([position_2x, position_2y, position_2z])
    position_3 = np.array([position_3x, position_3y, position_3z])
    
    # magnitudes of distances
    r1_2 = np.sqrt(np.sum((position_2 - position_1)**2))
    r1_3 = np.sqrt(np.sum((position_3 - position_1)**2))
    r2_3 = np.sqrt(np.sum((position_3 - position_2)**2))
    
    # right hand sides
    rhs = np.zeros(18)
    rhs[0] = velocity_1x
    rhs[1] = velocity_2x
    rhs[2] = velocity_3x
    rhs[3] = (G*M2*(position_2x - position_1x)/r1_2**3 + G*M3*(position_3x -
       position_1x)/r1_3**3)
    rhs[4] = (G*M3*(position_3x - position_2x)/r2_3**3 + G*M1*(position_1x -
       position_2x)/r1_2**3)
    rhs[5] = (G*M2*(position_2x - position_3x)/r2_3**3 + G*M1*(position_1x -
       position_3x)/r1_3**3)
    rhs[6] = velocity_1y
    rhs[7] = velocity_2y
    rhs[8] = velocity_3y
    rhs[9] = (G*M2*(position_2y - position_1y)/r1_2**3 + G*M3*(position_3y -
       position_1y)/r1_3**3)
    rhs[10] = (G*M3*(position_3y - position_2y)/r2_3**3 + G*M1*(position_1y -
       position_2y)/r1_2**3)
    rhs[11] = (G*M2*(position_2y - position_3y)/r2_3**3 + G*M1*(position_1y -
       position_3y)/r1_3**3)
    rhs[12] = velocity_1z
    rhs[13] = velocity_2z
    rhs[14] = velocity_3z
    rhs[15] = (G*M2*(position_2z - position_1z)/r1_2**3 + G*M3*(position_3z -
       position_1z)/r1_3**3)
    rhs[16] = (G*M3*(position_3z - position_2z)/r2_3**3 + G*M1*(position_1z -
       position_2z)/r1_2**3)
    rhs[17] = (G*M2*(position_2z - position_3z)/r2_3**3 + G*M1*(position_1z -
       position_3z)/r1_3**3)
    return rhs








N = 100
a = 0.0
b = 20.0
dt = 0.02
t_v = np.arange(0.0, 30.0, dt)
base_shape = (2,N)
base_qp_0 = np.zeros(base_shape, dtype=float)
base_qp_0[0,0] = np.pi/2.0







