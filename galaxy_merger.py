#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 10:05:02 2020

@author: blarsen10

Create model with no class:
    Controlled by single mass SMBH potential
    points in space represent stars
    create initial conditions that spin about star
    
Need to comment

"""

import matplotlib.animation as ani
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as spi

plt.close('all')

# =============================================================================
# Functions
# =============================================================================

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
    for i in range(num_galaxies):
        rhs[2*i] = state[2*i+1]
        rhs[2*i+1] = 0
    for i in range(N_total):
        rhs[2*i+2*num_galaxies] = state[2*i+2*num_galaxies+1] # velocity
        rhs[2*i+2*num_galaxies+1] = G*M*(state[0] - state[2*i+2*num_galaxies])/r_array[i]**3
    
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
    
    # calculate positions relative to bodies
    r_cent = np.sum((state[0] - state[2])**2)**0.5
    r1_array = np.zeros(N_total)
    r2_array = np.zeros(N_total)
    for i in range(N_total):
        r1_array[i] = np.sum((state[2*i+2*num_galaxies] - state[0])**2)**0.5
        r2_array[i] = np.sum((state[2*i+2*num_galaxies] - state[2])**2)**0.5
    
    # right hand sides (velocities and accelerations)
    rhs = np.zeros((num_equations//3, 3))
    rhs[0] = state[1]
    rhs[1] = G*M*(state[2] - state[0])/r_cent**3
    rhs[2] = state[3]
    rhs[3] = G*M*(state[0] - state[2])/r_cent**3
    for i in range(N_total):
        pos = state[2*i+2*num_galaxies]
        vel = state[2*i+2*num_galaxies+1]
        rhs[2*i+2*num_galaxies] = vel # velocity
        rhs[2*i+2*num_galaxies+1] = (G*M*(state[0] - pos)/r1_array[i]**3 +
                                     G*M*(state[2] - pos)/r2_array[i]**3)
    
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
    x_cent = np.zeros(num_galaxies)
    y_cent = np.zeros(num_galaxies)
    x = np.zeros(N_total)
    y = np.zeros(N_total)
    for j in range(num_galaxies):
        x_cent[j] = solution['y'][2*j,0][i]
        y_cent[j] = solution['y'][2*j,1][i]
    cent_bodies[0].set_data(x_cent, y_cent)
    for j in range(N_total):
        x[j] = solution['y'][2*j+num_galaxies*2, 0][i]
        y[j] = solution['y'][2*j+num_galaxies*2, 1][i]
    point_bodies[0].set_data(x, y)
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
num_galaxies = 2
galaxy_pos = np.array([[0, 0, 0], [5, 5, 0]])
galaxy_vel = np.array([[2, 0, 0], [-2, 0, 0]])

# rings of bodies with more on the edges
N_inner_ring = 10
N_outer_ring = 20
num_rings = 2
N_per_ring = np.linspace(N_outer_ring, N_inner_ring, num_rings, dtype = int)
N = np.sum(N_per_ring)
N_total = N*num_galaxies
r_outer = 2
check_N = True

if check_N:
    input("total N = " + str(N_total) + " press any enter to continue (control C to escape)")

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
t_max = 3 # years
Nt = 10000
times = np.linspace(t_min, t_max, Nt)

# other
animation_speed_scaling = 10
num_equations = ((N*num_galaxies)+num_galaxies)*6 # total number of equations tobe solved in solve_ivp

# =============================================================================
# Initial conditions
# =============================================================================

# initial conditions
init_cond = np.zeros((num_equations//3, 3))
for i in range(num_galaxies):
    init_cond[i*2] = galaxy_pos[i]
    init_cond[i*2+1] = galaxy_vel[i]
# build points in a circle (x and y are sin and cos)
index1 = num_galaxies*2
for q in range(num_galaxies):
    index2 = 0
    for i in range(num_rings):
        for j in range(N_per_ring[i]):
            init_cond[index1] = [r_array[index2]*np.sin(j/N_per_ring[i]*2*np.pi), r_array[index2]*np.cos(j/N_per_ring[i]*2*np.pi), 0] + init_cond[2*q]
            init_cond[index1+1] = [v_array[index2]*np.cos(j/N_per_ring[i]*2*np.pi), -v_array[index2]*np.sin(j/N_per_ring[i]*2*np.pi), 0] + init_cond[2*q+1]
            index1 += 2
            index2 += 1

# =============================================================================
# Solution
# =============================================================================

# solve_ivp yay
# it doesn't like a multideminsional init_cond array, so reshape it
# Radau to better conserve energy

func = {
    1: rhs_1_body,
    2: rhs_2_body
}
rhs = func[num_galaxies]
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
ani.FuncAnimation(fig, animate, init_func=init_animate, frames=1000,
                        interval=10, blit=True)
plt.xlim(-8.1, 13.1)
plt.ylim(-7.1, 11.1)

