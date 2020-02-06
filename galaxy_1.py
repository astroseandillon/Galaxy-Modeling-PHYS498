#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 10:05:02 2020

@author: blarsen10

Create model with no class:
    Controlled by single mass SMBH potential
    points in space represent stars
    create initial conditions that spin about star

"""

import matplotlib.animation as ani
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as spi

plt.close('all')

def rhs(time, state):
    
    # reshape state
    state = np.reshape(state, (4, 3))
    
    # single central mass, single particle surrounding
    # each state index contain x y and z data
    cent_pos = state[0]
    cent_vel = state[1]
    position_1 = state[2]
    velocity_1 = state[3]
    
    # magnitude of distance
    r1 = np.sqrt(np.sum((cent_pos - position_1)**2))
    
    # right hand sides (velocities and accelerations)
    rhs = np.zeros((4, 3))
    rhs[0] = cent_vel
    rhs[1] = 0
    rhs[2] = velocity_1
    rhs[3] = G*M*(cent_pos - position_1)/r1**3
    
    return np.reshape(rhs, 12)

# initialization function for animation
def init_animate():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

# what to animate at every frame
def animate(i):
    # scaling speeds up the animation
    scaling = 10
    i*=scaling
    # plot the trajectories at given frame
    x = solution['y'][2, 0][i]
    y = solution['y'][2, 1][i]
    line.set_data(x,y)
    # change the text to reflect current age
    time_text.set_text('years: ' + str(i*t_max/(1000*scaling)))
    return line, time_text

G = 4*np.pi**2 # AU^2/(year^2 * Msun)
M = 1 # Solar masses

t_min = 0 # years
t_max = 1 # years
N = 10000 # fixed 100 steps so animation isn't too slow
times = np.linspace(t_min, t_max, N)

# initial conditions
init_cond = np.zeros((4, 3))
init_cond[0] = [0, 0, 0]
init_cond[1] = [0, 0, 0]
init_cond[2] = [0, -1, 0]
# for circular orbit: v = (GMr)^0.5
init_cond[3] = [2*np.pi, 0, 0]

# solve_ivp yay
# it doesn't like a multideminsional init_cond array, so reshape it
# Radau to better conserve energy
solution = spi.solve_ivp(rhs, [t_min, t_max], np.reshape(init_cond, 12),
                          vectorized = False, t_eval = times, method='Radau')

# reshape again: x is state; y is x,y,z; z steps through time
solution['y'] = np.reshape(solution['y'], (4, 3, N))

fig = plt.figure()
# initialize points to represent bodies and text to be animated
line, = plt.plot([], [], '.k')
time_text = plt.text(0.15, 0.15, '', transform=plt.gcf().transFigure)
# central mass
plt.plot(solution['y'][0, 0, 0], solution['y'][0, 1, 0], 'ok')
# trajectory
#plt.plot(solution['y'][2, 0], solution['y'][2, 1], '--k')
plt.xlim(-1.1, 1.1)
plt.ylim(-1.1, 1.1)
# animate point masses
ani = ani.FuncAnimation(fig, animate, init_func=init_animate, frames=1000,
                  interval=10, blit = True)


