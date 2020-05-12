'''
Cindy Olvera Perez, Bjorn Larsen, Sean Dillon
April 22, 2020
Physics 498: Galaxy Simulation
'''

# In[0]
import sim, figure
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.animation as ani
from astropy.table import Table
from astropy.io import ascii

def begin(start,stop,timesteps,
        Vfilename, Sfilename, overwrite=True):
    '''
    starts the simulation and loops over according to the number of timesteps
    creates an astropy table and writes to a text file
    may 7, 2020: can only change the velocity in the j direction
    '''
#    j= np.arange(start,stop,timesteps)
    j=np.array([3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.8])
#    sol=[]
    for i in range(len(j)):
        if i==0:
            solution, num_galaxies, n_total, nt = sim.galaxy(gPos=np.array([[-2,-5,0],[2,5,0]]), 
                                                             gVel=np.array([[0,j[i],0],[0,-j[i],0]]))
            col= figure.label(len(j))
            ar= np.arange(0,2*len(j))
#            col_sol= figure.sLabels(nt+1)
            ar_even= ar[ar%2==0]
            ar_odd= ar[ar%2!=0]
            v, v_com, v_initial, v_final, v_cent, n_total = sim.velocities(num_galaxies, n_total, nt, solution)
            galaxyData= Table(data=[v_initial, v_final], names=[col[ar_even[i]], col[ar_odd[i]]])
#            sol.append(solution)
        else:
            solution, num_galaxies, n_total, nt = sim.galaxy(gVel=np.array([[0,j[i],0],[0,-j[i],0]]), 
                                                             gPos=np.array([[-2,-5,0],[2,5,0]]))
            v, v_com, v_initial, v_final, v_cent, n_total = sim.velocities(num_galaxies, n_total, nt, solution)
            figure.Vcolumn(galaxyData, v_initial, v_final, col[ar_even[i]], col[ar_odd[i]])
#            sol.append(solution)
#            figure.Scolumn(galaxySol, solution, col_sol[i])
#    col_sol= figure.sLabels(len(sol))
#    galaxySol = Table(sol, names=col_sol)
    ascii.write(galaxyData, Vfilename, format='commented_header', overwrite=overwrite)
#    ascii.write(galaxySol, Sfilename, format='commented_header', overwrite=overwrite)
#    return solution, num_galaxies


#start siulation
begin(start=4.8,stop=4.9,timesteps=0.2,
      Vfilename='galaxyModel3.txt', Sfilename='galaxySolutions1.txt')

# In[0]
#create histograms using textfile
figure.histograms('galaxyModel3.txt',binStart=5, binStop=17, vel= r'Velocities $\frac{AU}{year}$')


# In[1]
#figure.histograms('galaxyModel2.txt',binStart=5, binStop=17, vel= r'Velocities $\frac{AU}{year}$')
j = np.arange(4.8, 4.9, 0.2)
print(j)

