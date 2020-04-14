# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 11:29:20 2020

@author: Sean Dillon

"""
#from nelsons code
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')


def gen_force_pendulum(theta):
    #force is a function of position
    return -g/L*np.sin(theta)

def gen_momentum_pendulum(omega):
    #momentum is a function of velocity
    return omega


def simp_int_func(Force_func, Momentum_func, del_t, n_times, 
                  init_momentum, init_position):
    """
    This is the simplectic integrator function for all separable hamiltonian 
    equations.  
    
    Force_func:         Function that gives us force based on position
    Momentum_func:      Function that gives us momentum based on velocity
    del_t:              floating point number that tells us length of timestep
    n_times:            Int that describes total number of timesteps
    init_momentum:      initial momentum
    init_position:      initial position
    
    Returns:            Numpy array of type [position, momentum]
    """
    a1 = (2.0+2.0**(1.0/3.0)+2.0**(-1.0/3.0))/6.0
    a2 = (1.0 - 2.0**(1.0/3.0) - 2.0**(-1.0/3.0))/6.0
    a_coefs = np.array([a1, a2, a1, a2])
    b2 = 1.0/(2.0 - 2.0**(1.0/3.0))
    b3 = 1.0/(1.0 - 2.0**(2.0/3.0))
    b_coefs = np.array([0.0, b2, b3, b2])
    #integrator coefficients
    momentum = np.zeros(n_times)
    position = np.zeros(n_times)
    #creating our arrays
    position[0] = init_position
    momentum[0] = init_momentum
    #setting initial conditions
    for i in range(1,n_times):
        momentum_im1 = momentum[i-1]
        position_im1 = position[i-1]
        for j in range(4):
            momentum_i = momentum_im1 + b_coefs[j]*del_t*Force_func(position_im1)
            position_i = position_im1 + a_coefs[j]*del_t*Momentum_func(momentum_i)
            momentum_im1 = momentum_i
            position_im1 = position_i
        momentum[i] = momentum_i
        position[i] = position_i
    
    return np.array([position, momentum])

    

    
    
    
    
























