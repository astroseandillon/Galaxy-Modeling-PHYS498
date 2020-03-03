#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 10:05:02 2020
@author: blarsen10
Create model with no class:
    Controlled by single mass SMBH potential
    points in space represent stars
    create initial conditions that spin about star

March 3, 2020
Cindy Olvera 
Must run main.py before running animate.py
main function computes galaxy model & returns the following parameters
animation_speed_scaling, solution, t_max, N
animate.py takes the parameters above and creates an animation
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as spi


plt.close('all')

# =============================================================================
# Functions
# =============================================================================
def main(M=1, n=120, rings=5, outer_rings=1,      # constant member variables
           t_min=0, t_max=1                       # time
           
           ):
    def rhs(time, state):
        
        # reshape state
        state = np.reshape(state, (num_equations//3, 3))
        
        # single central mass, single particle surrounding
        # each state index contain x y and z data
        # cent_pos = state[0]
        # cent_vel = state[1]
        # position_1 = state[2]
        # velocity_1 = state[3]
        # position_2 = state[4]
        # and so on...
        
        # right hand sides (velocities and accelerations)
        rhs = np.zeros((num_equations//3, 3))
        rhs[0] = state[1]
        rhs[1] = 0
        for i in range(N):
            rhs[2*i+2] = state[2*i+3] # velocity
            rhs[2*i+3] = G*M*(state[0] - state[2*i+2])/r_array[i]**3
        
        return np.reshape(rhs, num_equations)
    
    # initialization function for animation
    def init_animate():
        line.set_data([], [])
        time_text.set_text('')
        return line, time_text
    
    # what to animate at every frame
    def animate(i):
        
        # scaling speeds up the animation
        i*=animation_speed_scaling
        
        # plot the trajectories at given frame
        x = np.zeros(N)
        y = np.zeros(N)
        for j in range(N):
            x[j] = solution['y'][2*j+2, 0][i]
            y[j] = solution['y'][2*j+2, 1][i]
        line.set_data(x,y)
        
        # change the text to reflect current age
        time_text.set_text('years: ' + str(i*t_max/(1000*animation_speed_scaling)))
        return line, time_text
    
    # =============================================================================
    # Parameters
    # =============================================================================
    
    '''central mass parameters'''
    G = 4*np.pi**2                               # AU^2/(year^2 * Mgal)
    M = M                                        # Solar masses
    
    # bodies; all at distance 1 AU away for now
    N = n
    
    ## rings of bodies with more on the edges
    #N_inner_ring = 12
    #N_outer_ring = 36
    num_rings = rings
    r_outer = outer_rings
    
    # get initial position and velocity values for each star
    r_array = np.zeros(N)
    N_per_ring = int(np.ceil(N/num_rings))
    for i in range(num_rings):
        first_index = i*N_per_ring
        last_index = np.min([N, (i+1)*N_per_ring])
        r_array[first_index:last_index] = r_outer*((num_rings-i)/num_rings)

    # velocity for circular orbit; positive for clockwise, negative for counterclockwise
    v_array = (G*M/r_array)**0.5
    
    # times
    t_min = t_min                                # years
    t_max = t_max                                # years
    Nt = 10000                                   # fixed 100 steps so animation isn't too slow
    times = np.linspace(t_min, t_max, Nt)
    
    # other
    animation_speed_scaling = 10
    num_equations = (N+1)*6                      # total number of equations tobe solved in solve_ivp
    
    # =============================================================================
    # Initial conditions
    # =============================================================================
    
    # initial conditions
    init_cond = np.zeros((num_equations//3, 3))
    init_cond[0] = [0, 0, 0]
    init_cond[1] = [0, 0, 0]
    
    # build points in a circle (x and y are sin and cos)
    for i in range(N):
        init_cond[2*i+2] = [r_array[i]*np.sin(i/N_per_ring*2*np.pi), r_array[i]*np.cos(i/N_per_ring*2*np.pi), 0]
        init_cond[2*i+3] = [v_array[i]*np.cos(i/N_per_ring*2*np.pi), -v_array[i]*np.sin(i/N_per_ring*2*np.pi), 0]
    
    # =============================================================================
    # Solution
    # =============================================================================
    
    '''
    solve_ivp yay
    it doesn't like a multideminsional init_cond array, so reshape it
    Radau to better conserve energy
    '''
    solution = spi.solve_ivp(rhs, [t_min, t_max], np.reshape(init_cond, num_equations),
                             t_eval = times, method='Radau')
    
    # reshape again: x is state; y is x,y,z; z steps through time
    solution['y'] = np.reshape(solution['y'], (num_equations//3, 3, Nt))
    
    return animation_speed_scaling, solution, t_max, N

animation_speed_scaling, solution, t_max, N = main()



