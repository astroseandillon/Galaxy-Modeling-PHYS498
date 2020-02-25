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
from integrator_exceptions import SalvagedResultException as exceptions

plt.close('all')



N=3
dt = 0.02
t_v = np.arange(0.0, 30.0, dt)
qp_0 = np.zeros((2,N), dtype=float)
qp_0[0,0] = np.pi/2.0





def V (q):
        """Potential energy is a function of the position only."""
        # np.linalg.norm(q) gives the angle from the vertical axis
        return -np.cos(np.linalg.norm(q))

def K (p):
        """Kinetic energy is a function of the momentum only.  It is assumed that the pendulum has unit mass."""
        return 0.5*np.sum(np.square(p))

def H (coordinates):
        """The Hamiltonian is the sum of kinetic and potential energy."""
        q = coordinates[0,:]
        p = coordinates[1,:]
        return K(p) + V(q)

def dK_dp (p):
        return p

def dV_dq (q):
    # sinc(x) is sin(pi*x)/(pi*x) when x is not 0 -- this is used to avoid the singularity in sin(x)/x.
    return -np.sinc(np.linalg.norm(q)/np.pi) * q

def dH_dq (q, p):
    return dV_dq(q)

def dH_dp (q, p):
    return p

UpdateStepCoefficients = collections.namedtuple('UpdateStepCoefficients', ['euler1', 'verlet2', 'ruth3', 'ruth4'])
update_step_coefficients = UpdateStepCoefficients(
    # euler1
    np.array([
        [1.0],
        [1.0]
    ]),
    # verlet2
    np.array([
        [0.0, 1.0],
        [0.5, 0.5]
    ]),
    # ruth3
    np.array([
        [1.0, -2.0/3.0, 2.0/3.0],
        [-1.0/24.0, 0.75, 7.0/24.0]
    ]),
    # ruth4
    __make_ruth4_update_step_coefficients()
)


a = integrate(
        initial_coordinates=qp_0,
        t_v=t_v,
        dK_dp=dK_dp(qp_0[1,:]),
        dV_dq = dV_dq(qp_0[0,:]),
        update_step_coefficients=update_step_coefficients.verlet2)














