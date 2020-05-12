#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 10:05:02 2020
@author: Bjorn, Cindy, Sean
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

#import csv
#import matplotlib.pyplot as plt
#import matplotlib as mpl
import numpy as np
#import os
import scipy.integrate as spi
import scipy.linalg as spl
import random
#from astropy.table import Column
#from astropy.io import ascii

# =============================================================================
# Functions
# =============================================================================
'''
Parameters:
    General galaxy and central mass parameters:
        g ; AU^2/(year^2 * Mgal)
        m ; Solar masses
        num_galaxies
        galaxy_pos ; initial position, array of shape (num_galaxies, 3)
        galaxy_vel ; initial velocity, array of shape (num_galaxies, 3)
        euler_angles ; Angles with which to rotate the galactic disc, array of
                       shape (num_galaxies, 3):
            1st val: angle with respect to x-axis (radians)
            2nd val: angle with respect to y-axis (radians)
            3rd val: angle with respect to z-axis (radians)
    Surrounding point body parameters:
        r_outer ; radius of galaxy/outermost ring
        n_inner ; number of bodies on innermost ring
        n_outer ; number of bodies on outermost ring
        num_rings
    Time data: (times in years)
        t_min
        t_max
        nt ; number of timesteps
    Other parameters:
        grav_gal ; set 0 for no mutual gravitation between galaxies (default 1)
        save ; saves solution and initial conditions as binary
        file_name ; name for the save data file
        dir_name ; name of folder to save data to
        check_n ; set true to output n before running
'''
def simulation(num_galaxies=1, 
         galaxy_pos=np.array([[0,0,0]]), galaxy_vel=np.array([[0,0,0]]), euler_angles=np.array([[0,0,0]]),
         r_outer=1, n_inner=12, n_outer=36, num_rings=5, 
         t_min=0, t_max=1,
         nt=1000, grav_gal=1, 
         check_n = True):
    # constant member variables 
    g = 4*np.pi**2
    m = 1
    
    def rhs_1_body(time, state):
        
        # reshape state
        state = np.reshape(state, (num_equations//3, 3))
        
        # single central mass, single particle surrounding
        # each state index contain x y and z data
        # cent_pos_1 = state[0]
        # cent_vel_1 = state[1]
        # cent_pos_2 = state[2]
        # etc...
        # position_1 = state[2*num_galaxies]
        # velocity_1 = state[2*num_galaxies+1]
        # position_2 = state[2*num_galaxies+2]
        # etc...
        
        # right hand sides (velocities and accelerations)
        rhs = np.zeros((num_equations//3, 3))
        for i in range(0, num_galaxies, 2):
            rhs[i] = state[i+1] # central velocity
            rhs[i+1] = 0 # central acceleration
        i_rhs = 2*num_galaxies # index for rhs and state arrays
        for i in range(n_total):
            rhs[i_rhs] = state[i_rhs+1] # velocities
            # a = GM(rhat)/r^2
            rhs[i_rhs+1] = (g*m*(state[0] - state[i_rhs])/r_array[i]**3)
            i_rhs += 2
        # reshape array back such that solve_ivp is happy with it
        return np.reshape(rhs, num_equations)
    
    def rhs_2_body(time, state):
        
        # reshape state
        state = np.reshape(state, (num_equations//3, 3))
        
        # single central mass, single particle surrounding
        # each state index contain x y and z data
        # cent_pos_1 = state[0]
        # cent_vel_1 = state[1]
        # cent_pos_2 = state[2]
        # etc...
        # position_1 = state[2*num_galaxies]
        # velocity_1 = state[2*num_galaxies+1]
        # position_2 = state[2*num_galaxies+2]
        # etc...
        
        # calculate positions of bodies to central masses
        r_cent = np.sum((state[0] - state[2])**2)**0.5
        r1_array = np.zeros(n_total)
        r2_array = np.zeros(n_total)
        # state index starts at 2*num_galaxies and iterates up by 2 each time
        i_state = 2*num_galaxies # index for rhs and state arrays
        for i in range(n_total):
            r1_array[i] = np.sum((state[i_state] - state[0])**2)**0.5
            r2_array[i] = np.sum((state[i_state] - state[2])**2)**0.5
            i_state += 2
        
        # right hand sides (velocities and accelerations)
        rhs = np.zeros((num_equations//3, 3))
        # central body 1
        rhs[0] = state[1]
        rhs[1] = grav_gal*g*m*(state[2] - state[0])/r_cent**3
        # central body 2
        rhs[2] = state[3]
        rhs[3] = grav_gal*g*m*(state[0] - state[2])/r_cent**3
        # massless bodies
        i_rhs = 2*num_galaxies # index for rhs and state arrays
        for i in range(n_total):
            pos = state[i_rhs]
            vel = state[i_rhs+1]
            rhs[i_rhs] = vel # velocity
            # a = GM(r1hat)/r1^2 + GM(r2hat)/r2^2
            if i_rhs < 4+n*2:
                rhs[i_rhs+1] = (g*m*(state[0] - pos)/r1_array[i]**3 + 
                                grav_gal*g*m*(state[2] - pos)/r2_array[i]**3)
            else:
                rhs[i_rhs+1] = (grav_gal*g*m*(state[0] - pos)/r1_array[i]**3 + 
                                g*m*(state[2] - pos)/r2_array[i]**3)
            i_rhs += 2
        
        # update where the solver is at
        PrintTime(time)
        
        return np.reshape(rhs, num_equations)

    '''
    MODIFIED:
    '''
    # This is another crude attempt to keep track of how the solver is progressing
    # ie using random numbers
    def PrintTime(time):
        if random.randint(1,500) == 1:
            print(str(np.round(time, 2)) + "/" + str(t_max) + " years calculated")

    # set up initial conditions
    def initial_conditions(r, v, angles):
        init_cond = np.zeros((num_equations//3, 3))
        # first set of initial conditions are for central bodies
        for i in range(num_galaxies):
            init_cond[i*2] = galaxy_pos[i]
            init_cond[i*2+1] = galaxy_vel[i]
        # define rotation matrices (1 per galaxy)
        # using angles alpha beta and gamma
        alpha = angles[:,0]
        beta = angles[:,1]
        gamma = angles[:,2]
        R = np.array([[np.cos(beta)*np.cos(gamma),
                       np.sin(alpha)*np.sin(beta)*np.cos(gamma)-np.cos(alpha)*np.sin(gamma),
                       np.cos(alpha)*np.sin(beta)*np.cos(gamma)+np.sin(alpha)*np.sin(gamma)],
                      [np.cos(beta)*np.sin(gamma),
                       np.sin(alpha)*np.sin(beta)*np.sin(gamma)+np.cos(alpha)*np.cos(gamma),
                       np.cos(alpha)*np.sin(beta)*np.sin(gamma)-np.sin(alpha)*np.cos(gamma)],
                      [-np.sin(beta),
                       np.sin(alpha)*np.cos(beta),
                       np.cos(alpha)*np.cos(beta)]])
        # build points in a circle (x and y are sin and cos)
        # it seems to require a lot of indexing
        i_ic = num_galaxies*2 # init cond index
        for gal_number in range(num_galaxies):
            i_rv = 0 # index of r and v arrays
            for i in range(num_rings):
                for j in range(n_per_ring[i]):
                    # build circular motion of points
                    init_cond[i_ic] = ([r[i_rv]*np.sin(j/n_per_ring[i]*2*np.pi),
                                     r[i_rv]*np.cos(j/n_per_ring[i]*2*np.pi), 0])
                    init_cond[i_ic+1] = ([-v[i_rv]*np.cos(j/n_per_ring[i]*2*np.pi),
                                       v[i_rv]*np.sin(j/n_per_ring[i]*2*np.pi), 0])
                    # multiply by rotation matrix to rotate points about axis
                    init_cond[i_ic] = spl.solve(R[:,:,gal_number], init_cond[i_ic])
                    init_cond[i_ic+1] = spl.solve(R[:,:,gal_number], init_cond[i_ic+1])
                    # add position/velocity of com
                    init_cond[i_ic] = init_cond[i_ic] + init_cond[2*gal_number]
                    init_cond[i_ic+1] = init_cond[i_ic+1] + init_cond[2*gal_number+1]
                    # continue weird indexing
                    i_ic += 2
                    i_rv += 1
        return init_cond

    # =============================================================================
    # Additional variables
    # =============================================================================

    # number of bodies on intermediate rings will be linear to the edges
    n_per_ring = np.linspace(n_outer, n_inner, num_rings, dtype = int)
    n = np.sum(n_per_ring) # calculate N
    n_total = n*num_galaxies # number of point bodies total
    
    # more time data
    times = np.linspace(t_min, t_max, nt)
    
    # total number of equations to be solved in solve_ivp
    num_equations = ((n*num_galaxies)+num_galaxies)*6
    
    # =============================================================================
    # Initial conditions
    # =============================================================================
    
    # get initial position and velocity values for each star
    if n > 0:
        # set up r values (around central mass) for each body
        r_array = np.zeros(n)
        r_vals = np.linspace(r_outer,r_outer/num_rings,num_rings)
        for i in range(num_rings):
            r_array[np.argmin(r_array):np.sum(n_per_ring[:i+1])] = r_vals[i]
        # velocity for circular orbit
        # positive for clockwise, negative for counterclockwise
        v_array = (g*m/r_array)**0.5
    
    # set up initial conditions
    init_cond = initial_conditions(r_array, v_array, euler_angles)

    # =============================================================================
    # Solution
    # =============================================================================
    
    # solve_ivp doesn't like vectorized arrays, so reshape it
    # Radau to better conserve energy
    func = {
        1: rhs_1_body,
        2: rhs_2_body
    }
    rhs = func[num_galaxies]
    print('Begin time evolution')
    solution = spi.solve_ivp(rhs, [t_min, t_max], np.reshape(init_cond, num_equations),
                             t_eval = times, method='Radau')
    
    # reshape again: x is state; y is x,y,z; z steps through time
    solution['y'] = np.reshape(solution['y'], (num_equations//3, 3, nt))
    
    return solution, num_galaxies, n_total, nt
'''
# returns velocities for point masses, point masses with respect to their central
# masses, and the velocities of the central masses
def velocities(ngal,n,nt,solution):
    v_cent = np.zeros((ngal,3,nt))
    v = np.zeros((n,3,nt))
    v_com = np.zeros((n,nt))
    for i in range(nt):
        # get central mass velocities
        for j in range(ngal):
            v_cent[j,:,i] = solution['y'][j*2+1,:,i]
        # get point mass velocities
        v[:,:,i] = solution['y'][ngal*2::2,:,i]
    # counter counts through stars in each galaxy
    count = 0
    gal = -1
    for i in range(n):
        # condition: counter has finished counting one galaxy
        if count % (n//ngal) == 0:
            gal += 1
        v_temp = v[i] - v_cent[gal]
        # rms velocities
        v_com[i] = np.sum(v_temp**2, axis=0)**0.5
        count += 1
    v_initial = v_com[:,0] # velocity before merger
    v_final = v_com[:,-1]   # velocity after merger 
    
    return v, v_com, v_initial, v_final, v_cent, n
'''
def velocities(ngal,n,nt,solution,extract_data=False,gal_num=0):
    v_cent = np.zeros((ngal,3,nt))
    v = np.zeros((n,3,nt))
    v_com = np.zeros((n,nt))
    for i in range(nt):
        # get central mass velocities
        for j in range(ngal):
            v_cent[j,:,i] = solution[j*2+1,:,i]
        # get point mass velocities
        v[:,:,i] = solution[ngal*2+1::2,:,i]
    # counter counts through stars in each galaxy
    count = 0
    gal = -1
    for i in range(n):
        # condition: counter has finished counting one galaxy
        if count % (n//ngal) == 0:
            gal += 1
        v_temp = v[i] - v_cent[gal]
        # rms velocities
        v_com[i] = np.sum(v_temp**2, axis=0)**0.5
        count += 1
    v_initial = v_com[:,0] # velocity before merger
    v_final = v_com[:,-1]   # velocity after merger 
#    if extract_data:
#        v = ExtractGalaxyData(v,ngal,gal_num)
#        v_com = ExtractGalaxyData(v_com,ngal,gal_num)
#        v_cent = v_cent[gal_num-1]
#        
    return v, v_com, v_cent

def galaxy(gVel, gPos):
    solution, num_galaxies, n_total, nt = simulation(num_galaxies=2, galaxy_pos=gPos,
                                   galaxy_vel=gVel,
                                   euler_angles=np.array([[0,0,0],[0,-np.pi/2,0]]),
                                   t_max=5, nt=2, r_outer=2, num_rings=1)
    return solution, num_galaxies, n_total, nt
