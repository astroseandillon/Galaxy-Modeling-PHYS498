'''
Cindy Olvera Perez, Bjorn, Sean Dillon
April 22, 2020
Physics 498: Galaxy Simulation
'''

# In[0]
import main1 as main
import numpy as np 
from astropy.table import Table


# galaxy merger with galactic discs in the yz plane
solution, num_galaxies, n_total, nt = main.main(num_galaxies=2, galaxy_pos=np.array([[-2,-5,0],[2,5,0]]),
                                   galaxy_vel=np.array([[0,3,0],[0,-3,0]]),
                                   euler_angles=np.array([[0,0,0],[0,np.pi/2,0]]),
                                   t_max=5, nt=2001, r_outer=2, num_rings=1)

# In[1]

v, v_com, v_initial, v_final, v_cent, n_total = main.velocities(num_galaxies, n_total, nt, solution)
title1 = ""
xlabel = 'Velocities (AU/year)'

main.galaxyHistogram(v_initial, v_final, n_total, title1, xlabel, binw=0.2)


# In[2]
# saving data
'''
index = np.arange(0, n_total)

galaxyData = Table([index], names='i')

main.column(galaxyData, v_initial, "initial velocity 1")

print(galaxyData)
'''
