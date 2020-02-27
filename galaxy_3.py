#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 10:05:02 2020

@author: blarsen10

TO DO:
    CLEAN UP CODE
    Classes/Numba/Decorators for function time-output
    *Fix animation saving
    Symplectic integrator
    Allow for unequal galaxies
    Adaptive timesteps?
    Try other initial conditions
    Make it a true N-body code
    Energy analysis
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
    for i in range(0, num_galaxies, 2):
        rhs[i] = state[i+1] # central velocity
        rhs[i+1] = 0 # central acceleration
    i_rhs = 2*num_galaxies # index for rhs and state arrays
    for i in range(N_total):
        rhs[i_rhs] = state[i_rhs+1] # velocities
        # a = GM(rhat)/r^2
        rhs[i_rhs+1] = (G*M*(state[0] - state[i_rhs])/r_array[i]**3)
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
    r1_array = np.zeros(N_total)
    r2_array = np.zeros(N_total)
    # state index starts at 2*num_galaxies and iterates up by 2 each time
    i_state = 2*num_galaxies # index for rhs and state arrays
    for i in range(N_total):
        r1_array[i] = np.sum((state[i_state] - state[0])**2)**0.5
        r2_array[i] = np.sum((state[i_state] - state[2])**2)**0.5
        i_state += 2
    
    # right hand sides (velocities and accelerations)
    rhs = np.zeros((num_equations//3, 3))
    # central body 1
    rhs[0] = state[1]
    rhs[1] = G*M*(state[2] - state[0])/r_cent**3
    # central body 2
    rhs[2] = state[3]
    rhs[3] = G*M*(state[0] - state[2])/r_cent**3
    # massless bodies
    i_rhs = 2*num_galaxies # index for rhs and state arrays
    for i in range(N_total):
        pos = state[i_rhs]
        vel = state[i_rhs+1]
        rhs[i_rhs] = vel # velocity
        # a = GM(r1hat)/r1^2 + GM(r2hat)/r2^2
        rhs[i_rhs+1] = (G*M*(state[0] - pos)/r1_array[i]**3 + 
                        G*M*(state[2] - pos)/r2_array[i]**3)
        i_rhs += 2
    
    # update where the solver is at
    PrintTime(time)
    
    return np.reshape(rhs, num_equations)

# This is a crude attempt to keep track of how the solver is progressing
def PrintTime(time):
    global t_count
    if t_count % Nt == 0:
        print(str(np.round(time, 2)) + "/" + str(t_max) + " years calculated")
    t_count += 1

# set up initial conditions
def initial_conditions():
    init_cond = np.zeros((num_equations//3, 3))
    # first set of initial conditions are for central bodies
    for i in range(num_galaxies):
        init_cond[i*2] = galaxy_pos[i]
        init_cond[i*2+1] = galaxy_vel[i]
    # build points in a circle (x and y are sin and cos)
    # it seems to require a lot of indexing
    i_ic = num_galaxies*2 # init cond index
    for gal_number in range(num_galaxies):
        i_rv = 0 # index of r and v arrays
        for i in range(num_rings):
            for j in range(N_per_ring[i]):
                init_cond[i_ic] = ([r_array[i_rv]*np.sin(j/N_per_ring[i]*2*np.pi),
                                 r_array[i_rv]*np.cos(j/N_per_ring[i]*2*np.pi),
                                 0] + init_cond[2*gal_number])
                init_cond[i_ic+1] = ([v_array[i_rv]*np.cos(j/N_per_ring[i]*2*np.pi),
                                   -v_array[i_rv]*np.sin(j/N_per_ring[i]*2*np.pi),
                                   0] + init_cond[2*gal_number+1])
                i_ic += 2
                i_rv += 1
    return init_cond

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
    # set up data pairs for animations of all point bodies
    j_sol = num_galaxies*2
    for j in range(N_total):
        x[j] = solution['y'][j_sol, 0][i]
        y[j] = solution['y'][j_sol, 1][i]
        j_sol += 2
    # loop across all galaxies
    last_index = 0
    for j in range(num_galaxies):
        x_cent[j] = solution['y'][2*j,0][i]
        y_cent[j] = solution['y'][2*j,1][i]
        # update the artists; the complicated indexing is just to identify
        # which point bodies are associated with which galaxy
        cent_bodies[j].set_data(x_cent[j], y_cent[j])
        next_index = N_total*(j+1)//num_galaxies
        point_bodies[j].set_data(x[last_index:next_index],
                                 y[last_index:next_index])
        last_index = next_index
    # change the text to reflect current age
    time_text[0].set_text('years: ' + str(i*t_max/Nt))
    return patches

# =============================================================================
# Parameters
# =============================================================================

# general parameters
G = 4*np.pi**2 # AU^2/(year^2 * Mgal)
M = 1 # Solar masses

# SET UP FOR 2 EQUAL GALAXIES
num_galaxies = 2
# galaxy intial position and velocity
galaxy_pos = np.array([[0, 0, 0], [5, 5, 0]])
galaxy_vel = np.array([[1.5, 0, 0], [-1.5, 0, 0]])

# rings of bodies with more on the edges
r_outer = 1 # radius of galaxy
N_inner_ring = 5 # number of bodies on innermost ring
N_outer_ring = 10 # number of bodies on outermost ring
num_rings = 2
# number of bodies on intermediate rings will be linear to the edges
N_per_ring = np.linspace(N_outer_ring, N_inner_ring, num_rings, dtype = int)
N = np.sum(N_per_ring) # calculate N
N_total = N*num_galaxies # number of point bodies total

# times
t_min = 0 # years
t_max = 1 # years
Nt = 1000 # number of timesteps (consider defining dt?)
t_interval = Nt//t_max # number of timesteps in one year
t_count = 0 # Used in outputting time data during runtime
times = np.linspace(t_min, t_max, Nt)

# other parameters
check_N = True # output number of Ns before running
animation_speed_scaling = 10 # increase for faster animation
save_animation = False # animation saving is still broken
animation_writer = 'imagemagick' # use either this or ffmpeg, whichever is installed
animation_dir = "animations" # save animations to this directory
# if create custom name, do not start with 'ani'
animation_file_name = '' # leave blank to automatically create file name
# total number of equations tobe solved in solve_ivp
num_equations = ((N*num_galaxies)+num_galaxies)*6

# =============================================================================
# Initial conditions setup
# =============================================================================

# confirm you are not breaking your computer before running further
if check_N:
    input("total N = " + str(N_total) +
          " press any enter to continue (control C to escape)")

# get initial position and velocity values for each star
if N > 0:
    # set up r values (around central mass) for each body
    r_array = np.zeros(N)
    r_vals = np.linspace(r_outer,r_outer/num_rings,num_rings)
    for i in range(num_rings):
        r_array[np.argmin(r_array):np.sum(N_per_ring[:i+1])] = r_vals[i]
    # velocity for circular orbit
    # positive for clockwise, negative for counterclockwise
    v_array = (G*M/r_array)**0.5

# set up initial conditions
init_cond = initial_conditions()

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
solution = spi.solve_ivp(rhs, [t_min, t_max], np.reshape(init_cond, num_equations),
                         t_eval = times, method='Radau')

# reshape again: x is state; y is x,y,z; z steps through time
solution['y'] = np.reshape(solution['y'], (num_equations//3, 3, Nt))

# =============================================================================
# Plotting
# =============================================================================

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
animation = ani.FuncAnimation(fig, animate, init_func=init_animate,
                              frames=Nt//animation_speed_scaling,
                              interval=10, blit=True)
# Save animation
if save_animation:
    # create new animation file name, if one does not exist already
    if animation_file_name == '':
        files = ([f for f in os.listdir(animation_dir)
                 if os.path.isfile(os.path.join(animation_dir, f))])
        file_num = np.zeros(len(files))
        for file in files:
            if file[:3] == 'ani':
                extension = int(file[3:-4])
        if len(files) == 0:
            animation_file_name = 'ani1.mp4'
        else:
            animation_file_name = 'ani' + str(np.max(extension) + 1) + '.mp4'
    # save animation using imagemagick
    # currently horribly broken but it's trying its best
    animation.save(animation_dir + '/' + animation_file_name,
                   writer=animation_writer, fps=15, bitrate=1800,
                   metadata=dict(artist='Chico_Astro'))
plt.xlim(-8.1, 13.1)
plt.ylim(-7.1, 11.1)

