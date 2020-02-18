#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 10:05:02 2020

@author: blarsen10

Create model with no class:
    Controlled by single mass SMBH potential
    points in space represent stars
    create initial conditions that spin about star
    
Need to determine if difficulties conserving momentum are my fault of solve_ivp's fault

"""

import matplotlib.animation as ani
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as spi

plt.close('all')

# =============================================================================
# Functions
# =============================================================================

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
    cent_bodies[0].set_data([], [])
    point_bodies[0].set_data([], [])
    time_text[0].set_text('')
    return patches

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
    cent_bodies[0].set_data(solution['y'][0,0][i], solution['y'][0,1][i])
    point_bodies[0].set_data(x,y)
    # change the text to reflect current age
    time_text[0].set_text('years: ' + str(i*t_max/(1000*animation_speed_scaling)))
    return patches

# =============================================================================
# Parameters
# =============================================================================

# central mass parameters
G = 4*np.pi**2 # AU^2/(year^2 * Mgal)
M = 1 # Solar masses

# SET UP FOR 2 EQUAL GALAXIES
num_galaxies = 1

# rings of bodies with more on the edges
N_inner_ring = 12
N_outer_ring = 36
num_rings = 5
N_per_ring = np.linspace(N_outer_ring, N_inner_ring, num_rings, dtype = int)
N = np.sum(N_per_ring)
r_outer = 1
check_N = False

if check_N:
    input("total N = " + str(N) + " press any enter to continue (control C to escape)")

# get initial position and velocity values for each star
if N > 0:
    r_array = np.zeros(N)
    r_vals = np.linspace(r_outer,r_outer/num_rings,num_rings)
    for i in range(num_rings):
        r_array[np.argmin(r_array):np.sum(N_per_ring[:i+1])] = r_vals[i]
    # velocity for circular orbit; positive for clockwise, negative for counterclockwise
    v_array = (G*M/r_array)**0.5

# times
t_min = 0 # years
t_max = 1 # years
Nt = 10000 # fixed 100 steps so animation isn't too slow
times = np.linspace(t_min, t_max, Nt)

# other
animation_speed_scaling = 10
num_equations = (N+1)*6 # total number of equations tobe solved in solve_ivp

# =============================================================================
# Initial conditions
# =============================================================================

# initial conditions
init_cond = np.zeros((num_equations//3, 3))
init_cond[0] = [0, 0, 0]
init_cond[1] = [1, 0, 0]
# build points in a circle (x and y are sin and cos)
index = 0
for i in range(num_rings):
    for j in range(N_per_ring[i]):
        init_cond[2*index+2] = [r_array[index]*np.sin(j/N_per_ring[i]*2*np.pi), r_array[index]*np.cos(j/N_per_ring[i]*2*np.pi), 0] + init_cond[0]
        init_cond[2*index+3] = [v_array[index]*np.cos(j/N_per_ring[i]*2*np.pi), -v_array[index]*np.sin(j/N_per_ring[i]*2*np.pi), 0] + init_cond[1]
        index += 1

# =============================================================================
# Solution
# =============================================================================

# solve_ivp yay
# it doesn't like a multideminsional init_cond array, so reshape it
# Radau to better conserve energy
solution = spi.solve_ivp(rhs, [t_min, t_max], np.reshape(init_cond, num_equations),
                         t_eval = times, method='Radau')

# reshape again: x is state; y is x,y,z; z steps through time
solution['y'] = np.reshape(solution['y'], (num_equations//3, 3, Nt))

# =============================================================================
# Plotting
# =============================================================================

fig = plt.figure()
# initialize points to represent bodies and text to be animated
cent_bodies, = [plt.plot([], [], 'ok')] # central bodies
point_bodies, = [plt.plot([], [], '.k')] # point bodies
time_text = [plt.text(0.15, 0.15, '', transform=plt.gcf().transFigure)]
patches = cent_bodies + point_bodies + time_text
# animate
ani = ani.FuncAnimation(fig, animate, init_func=init_animate, frames=1000,
                        interval=10, blit=True)
plt.xlim(-1.1, 2.1)
plt.ylim(-1.1, 1.1)

