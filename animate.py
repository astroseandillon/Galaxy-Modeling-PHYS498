# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 10:28:49 2020

@author: Cindy
"""

import matplotlib.animation as ani
import matplotlib.pyplot as plt
import numpy as np

plt.close('all')

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

fig = plt.figure()
# initialize points to represent bodies and text to be animated
line, = plt.plot([], [], '.k')
time_text = plt.text(0.15, 0.15, '', transform=plt.gcf().transFigure)
# central mass
plt.plot(solution['y'][0, 0, 0], solution['y'][0, 1, 0], 'ok')
# animate point masses
ani = ani.FuncAnimation(fig, animate, init_func=init_animate, frames=1000,
                  interval=10, blit = True)
plt.xlim(-1.1, 1.1)
plt.ylim(-1.1, 1.1)