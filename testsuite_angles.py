# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 01:51:34 2020

@author: Bjorn

Who knows if this will actually work
"""

import main5 as main
import numpy as np
import random

# change rotation: z-axis to y-axis to x-axis to z-axis.
dir_name = 'data/angles_data_1/regular'
num_angle_changes = 12
max_angle = np.pi/2
angle_change = max_angle/num_angle_changes
angle_range = np.linspace(0, max_angle - angle_change, num_angle_changes)
angles = np.zeros((num_angle_changes*3, 3))
# build angles
angles[:12,0] = angle_range
angles[12:24,2] = angle_range
angles[12:24,0] += np.pi/2
angles[24:,1] = angle_range
angles[24:,0] += np.pi/2
angles[24:,2] += np.pi/2

for i in range(9,len(angles)):
    file_name = 'reg_run_' + str(i+1)
    print('starting for ' + file_name)
    main.main(num_galaxies=2, galaxy_pos=np.array([[-2,-5,0],[2,5,0]]),
              galaxy_vel=np.array([[0,3,0],[0,-3,0]]), 
              euler_angles=np.array([[0,0,0],angles[i]]), r_outer=2,
              t_max=5, nt=2001, dir_name=dir_name, file_name=file_name,
              check_n=False)

# change rotation: Monte Carlo style
dir_name = 'data/angle_data_1/mc'
num_trials = 20

# get random angles
angles = np.zeros(num_trials*3)
for i in range(num_trials*3):
    angles[i] = random.uniform(-np.pi/2,np.pi/2)
angles = np.reshape(angles, (num_trials,3))

for i in range(len(angles)):
    file_name = 'rand_run_' + str(i)
    print('starting for ' + file_name)
    main.main(num_galaxies=2, galaxy_pos=np.array([[-2,-5,0],[2,5,0]]),
              galaxy_vel=np.array([[0,3,0],[0,-3,0]]), 
              euler_angles=np.array([[0,0,0],angles[i]]), r_outer=2,
              t_max=5, nt=2001, dir_name=dir_name, file_name=file_name,
              check_n=False)



