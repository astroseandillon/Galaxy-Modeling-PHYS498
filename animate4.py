# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 10:28:49 2020

@author: Cindy, Bjorn, Sean

- Need to figure out why rhs_2_body takes so long to run every time
    - investigate by running on my computer vs the 123 computer
- Need to adapt the right hand side function in main to work with the integrator function
- Need to create option to load solution from data file and animate rather than running
the simulation again
"""

import main4 as main # rename as needed
import matplotlib.animation as ani
import matplotlib.pyplot as plt
import numpy as np

plt.close('all')

'''
    Animation parameters:
        x_spacing ; additional space to leave on plots in x-direction
        y_spacing ; additional space to leave on plots in y-direction
        animation_speed_scaling ; increase for faster animation
        save_animation
        animation_writer ; recommended either imagemagick or ffmpeg
        animation_dir ; save animations to this directory
        animation_file_name ; option to create custom file name:
            - .gif for imagemagick, .mp4 for ffmpeg
            - leave blank to automatically create file name
            - do not start with '_ani' unless you want to confuse my dumb program
'''

x_spacing = 0.1
y_spacing = 0.1
animation_speed_scaling = 1
save_animation = False
animation_writer = 'imagemagick'
animation_dir = 'animations' # make this folder if it doesn't exist already
animation_file_name = ''
histograms = False

# I am having the animate function call main so all work done toying with
# parameters and initial conditions can be done just in this file

# some fun examples below
# rotated galaxy; euler angles in radians
solution, num_galaxies = main.main(num_galaxies=1, galaxy_pos=np.array([[-10,-2,0]]),
                                   galaxy_vel=np.array([[0,0,0]]), euler_angles=np.array([[1.47,0.1,0.3]]),
                                   t_max=1.0, nt=1001, num_rings=5, check_n=True)
'''
# multiple galaxies at different rotations with no gravity
# also not sure why this is taking long time to run? not much is going on here
solution, num_galaxies = main.main(num_galaxies=2, galaxy_pos=np.array([[-1,0,0],[1,0,0]]),
                                   galaxy_vel=np.array([[0,0,0],[0,0,0]]),
                                   euler_angles=np.array([[1.47,0.1,0.3],[np.pi/2,-0.3,0]]),
                                   t_max=0.1, nt=201, num_rings=3, check_n=True, grav_gal=0)
# galaxy merger with galactic discs in the yz plane
solution, num_galaxies = main.main(num_galaxies=2, galaxy_pos=np.array([[-2,-5,0],[2,5,0]]),
                                   galaxy_vel=np.array([[0,3,0],[0,-3,0]]),
                                   euler_angles=np.array([[0,np.pi/2,0],[0,np.pi/2,0]]),
                                   t_max=5, nt=2001, r_outer=2, num_rings=1,
                                   check_n=True)
'''




n_total = np.shape(solution['y'])[0]//2-num_galaxies
t_max = np.max(solution['t'])
nt = np.size(solution['t'])

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
    x = np.zeros(n_total)
    y = np.zeros(n_total)
    # set up data pairs for animations of all point bodies
    j_sol = num_galaxies*2
    for j in range(n_total):
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
        next_index = n_total*(j+1)//num_galaxies
        point_bodies[j].set_data(x[last_index:next_index],
                                 y[last_index:next_index])
        last_index = next_index
    # change the text to reflect current age
    time_text[0].set_text('years: ' + format(solution['t'][i],'.2'))
    return patches

def plot_bounds():
    x_min = np.min(solution['y'][::2,0,:]) - x_spacing
    x_max = np.max(solution['y'][::2,0,:]) + x_spacing
    y_min = np.min(solution['y'][::2,1,:]) - y_spacing
    y_max = np.max(solution['y'][::2,1,:]) + y_spacing
    return np.array([x_min, x_max, y_min, y_max])

def SaveAnimation(file_name, writer, fps=60, bitrate=-1):
    print('Saving animation')
    # create new animation file name, if one does not exist already
    if file_name == '':
        file_name = main.FileName('ani', animation_dir, data_type = '.gif')
    # save animation as gif using imagemagick
    animation.save(animation_dir + '/' + file_name,
                   writer=writer, fps=fps, bitrate=bitrate,
                   metadata=dict(artist='Chico_Astro'))

'''
PLOTTING
'''
plt.close('all')
# galaxy animation
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
                              frames=nt//animation_speed_scaling,
                              interval=10, blit=True)
# determine plot bounds
bounds = plot_bounds()
plt.xlim(bounds[0], bounds[1])
plt.ylim(bounds[2], bounds[3])
plt.xlabel('x (AU)')
plt.ylabel('y (AU)')
# Save animation
if save_animation:
    SaveAnimation(animation_file_name, animation_writer, fps=60, bitrate=-1)
plt.show()

'''
WIP
'''
if histograms:
    # histogram of velocity distributions with respect to center of mass for each galaxy
    v_cent, v, v_com = main.velocities(num_galaxies, n_total, nt, solution)
    
    num_bins = n_total//4
    plt.figure()
    plt.title("velocity distribution before merger")
    plt.xlabel('velocities (AU/year)')
    plt.hist(v_com[:,0],num_bins)
    
    plt.figure()
    plt.title("velocity distribution after merger")
    plt.xlabel('velocities (AU/year)')
    plt.hist(v_com[:,-1],num_bins)
