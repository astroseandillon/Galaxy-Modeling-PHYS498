#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 10:05:02 2020

@author: blarsen10

TO DO:
    CLEAN UP CODE
    Symplectic integrator
    Adaptive timesteps?
    Try other initial conditions
    Make it a true N-body code
"""

import matplotlib.animation as ani
import matplotlib.pyplot as plt
import numpy as np
import os
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
    
    PrintTime(time)
    
    return np.reshape(rhs, num_equations)

def PrintTime(time):
    global t_count
    if t_count % t_interval == 0:
        print(str(int(time)) + "/" + str(t_max) + " years calculated")
    t_count += 1
'''
t_count = 0
t_min = 0
t_max = 15
Nt = 30000
t_interval = Nt//t_max
time = np.linspace(t_min, t_max, Nt)
for i in range(Nt):
    PrintTime(time[i])'''

# initialization function for animation
def init_animate():
    for i in range(num_galaxies):
        cent_bodies[i].set_data([], [])
        point_bodies[i].set_data([], [])
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
    for j in range(N_total):
        x[j] = solution['y'][2*j+num_galaxies*2, 0][i]
        y[j] = solution['y'][2*j+num_galaxies*2, 1][i]
    last_index = 0
    for j in range(num_galaxies):
        x_cent[j] = solution['y'][2*j,0][i]
        y_cent[j] = solution['y'][2*j,1][i]
        # update the artists
        cent_bodies[j].set_data(x_cent[j], y_cent[j])
        next_index = N_total*(j+1)//num_galaxies
        point_bodies[j].set_data(x[last_index:next_index], y[last_index:next_index])
        last_index = next_index
    # change the text to reflect current age
    time_text[0].set_text('years: ' + str(i*t_max/Nt))
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
galaxy_vel = np.array([[1.5, 0, 0], [-1.5, 0, 0]])

# rings of bodies with more on the edges
N_inner_ring = 5
N_outer_ring = 15
num_rings = 3
N_per_ring = np.linspace(N_outer_ring, N_inner_ring, num_rings, dtype = int)
N = np.sum(N_per_ring)
N_total = N*num_galaxies
r_outer = 1.5
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
t_max = 15 # years
Nt = 30000
t_interval = Nt//t_max
t_count = 0 # for use in outputting time data during solve_ivp runtime
times = np.linspace(t_min, t_max, Nt)

# other
animation_speed_scaling = 10
animation_dir = "animations"
save_animation = False # do not check unless ffmpeg is installed
file_name = '' # leave blank to automatically create file name; if create custom name, do not start with 'ani'
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

if save_animation:
    Writer = ani.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Chico_Astro'), bitrate=1800)
    if file_name != '':
        files = [f for f in os.listdir(animation_dir) if os.path.isfile(os.path.join(animation_dir, f))]
        file_num = np.zeros(len(files))
        for file in files:
            if file[:3] == 'ani':
                extension = int(file[3:-4])
        if len(files) == 0:
            file_name = 'ani1.mp4'
        else:
            file_name = 'ani' + str(np.max(extension) + 1) + '.mp4'

fig = plt.figure()
# initialize points to represent bodies and text to be animated
# central bodies
cent_body1, = [plt.plot([], [], 'ok')]
cent_body2, = [plt.plot([], [], 'or')]
cent_bodies = cent_body1 + cent_body2
# point bodies
point_bodies1, = [plt.plot([], [], '.k')]
point_bodies2, = [plt.plot([], [], '.r')]
point_bodies = point_bodies1 + point_bodies2
# text
time_text = [plt.text(0.15, 0.15, '', transform=plt.gcf().transFigure)]
patches = cent_bodies + point_bodies + time_text
# animate
animation = ani.FuncAnimation(fig, animate, init_func=init_animate, frames=Nt//animation_speed_scaling,
                        interval=10, blit=True)
if save_animation:
    animation.save(animation_dir + '/' + file_name)
plt.xlim(-8.1, 13.1)
plt.ylim(-7.1, 11.1)

