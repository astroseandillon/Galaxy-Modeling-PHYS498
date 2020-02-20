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
from three_body_problem import rhs
import matplotlib.pyplot as plt



def f(x):
    return x**4+3*x**4 + 5

N = 100
a = 0.0
b = 20.0
dt = 0.02
t_v = np.arange(0.0, 30.0, dt)
base_shape = (2,N)
base_qp_0 = np.zeros(base_shape, dtype=float)
base_qp_0[0,0] = np.pi/2.0







