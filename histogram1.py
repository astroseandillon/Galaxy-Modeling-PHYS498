'''
Cindy Olvera Perez, Bjorn, Sean Dillon
April 22, 2020
Physics 498: Galaxy Simulation
'''

# In[0]
import importlib
main = importlib.import_module('main5')
import numpy as np 
import matplotlib.pyplot as plt
from astropy.table import Table

# In[1]

num_runs = 36
directory = "regular"
run = "reg_run"


file_name = "data_files/angles_data/"+directory+"/"+run+"_" + str(i)
file = open(file_name,'rb')
data = np.load(file,allow_pickle=True)
solution = data['solution']
parameters = data['parameters']
file.close()


# In[2]

plt.close('all')

g = parameters[0]
m = parameters[1]
num_galaxies = parameters[2]
angles = np.zeros(3)
for j in range(3):
    angles[j] = round(parameters[5][1][j],2)
n_total = int(np.size(solution, axis=0)/2 - num_galaxies)
nt = np.size(solution, axis=2)

r, r_com, r_cent = main.positions(num_galaxies, n_total, nt, solution,
                                  extract_data=True, gal_num=2)
v, v_com, v_cent = main.velocities(num_galaxies, n_total, nt, solution,
                                   extract_data=True, gal_num=2)
n = len(r)

e_div_m = 0.5*v_com**2 - g*m/r_com


# Positions plot
title1 = "Galaxy 2 Relative Positions (angles: " + str(angles) +")"
label_i = "Position initial"
label_f = "Position final"
xlabel = 'Positions (AU)'
ylabel = 'Number Density ($\\frac{n(r)}{N_{total}}$)'
r_initial = r_com[:,0]
r_final = r_com[:,-1]
save_name = 'figures/'+directory+'/pos_' + str(i)
main.galaxyHistogram(r_initial, r_final, n, title1, xlabel, ylabel,
                     label_i, label_f, save=True, save_name=save_name)

# Velocities plot
title1 = "Galaxy 2 Velocities (angles: " + str(angles) +")"
label_i = "Velocity initial"
label_f = "Velocity final"
xlabel = 'Velocities ($\\frac{AU}{year}$)'
ylabel = 'Number Density ($\\frac{n(v)}{N_{total}}$)'
v_initial = v_com[:,0]
v_final = v_com[:,-1]
save_name = 'figures/'+directory+'/vel_' + str(i)
main.galaxyHistogram(v_initial, v_final, n, title1, xlabel, ylabel,
                     label_i, label_f, save=True, save_name=save_name)

# Energies plot
title1 = "Galaxy 2 Total Energies per unit mass (angles: " + str(angles) +")"
label_i = "Energy initial"
label_f = "Energy final"
xlabel = 'Energies per unit mass ($\\frac{AU}{year})^2$'
ylabel = 'Number Density ($\\frac{n(v)}{N_{total}}$)'
e_initial = e_div_m[:,0]
e_final = e_div_m[:,-1]
save_name = 'figures/'+directory+'/en_' + str(i)
main.galaxyHistogram(e_initial, e_final, n, title1, xlabel, ylabel,
                     label_i, label_f, binw=1, save=True, save_name=save_name)


