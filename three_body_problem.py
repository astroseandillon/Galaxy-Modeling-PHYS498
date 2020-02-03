# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

Author: Bjorn Larsen

Now that you’ve done two gravitating bodies in two dimensions, let’s have a little more
fun and look at the interaction of three gravitating bodies in three dimensions. This is a
seemingly simple problem, but it turns out to have some extremely rich dynamics up to
and including dynamical chaos. Let’s assume you have three stars which all have the mass
of the Sun and which only interact via gravity. The equations of motion for these stars are
given by

    d2x1/dt2 = G*M2*(x2 - x1)/|x2 - x1|^3 + G*M3*(x3 - x1)/|x3 - x1|^3
    d2x2/dt2 = G*M3*(x3 - x2)/|x3 - x2|^3 + G*M1*(x1 - x2)/|x1 - x2|^3
    d2x3/dt2 = G*M2*(x2 - x3)/|x2 - x3|^3 + G*M1*(x1 - x3)/|x1 - x3|^3

where ~x1, ~x2, ~x3 are the vector positions of stars 1, 2, and 3 as functions of time. As
was the case previously, we will need initial positions and velocities, for a total of eighteen
initial conditions.

These equations should conserve energy and momentum. The total energy of the system
is given by

    Etot = (1/2)*(M1*v1^2 + M2*v2^2 + M3*v3^2) - G*(M1*M2/|x1 - x2| + M1*M3/|x1 - x3| + M2*M3/|x2 - x3|)

The total momentum of the system should be given by

    ptot = M1*dx1/dt + M2*dx2/dt + M3*dx3/dt

As you learned in introductory mechanics, the total energy and all three components of
the total momentum should be conserved over time.
Construct evolution equations for each component of the position vectors and then
make them all first order (eighteen equations total). Be sure to include these equations
in the comments at the top of your code. Use scipy.integrate.solve ivp in the SciPy
library to evolve these equations in time. Explore the evolution for different combinations

of masses and initial conditions until you find a system where all three bodies start within
1 AU of each other and remain bound for at least 100 years. Develop a useful way to
display the results. Also plot the total energy and momentum of your system as a function of time.
What do you notice? Finally, explore several different methods available in
solve ivp using the method keyword.

Pseudocode:
    imports
    functions:
        right hand side function:
            6 vector equations expanded into 18 when I realized solve_ivp does
            not like vector equations
            6 state equations:
                position_1 (array of size 3)
                position_2
                position_3 
                velocity_1
                velocity_2
                velocity_3
            6 first order differential equations representing two coupled second order equations:
                rhs_1 = velocity_1
                rhs_2 = velocity_2
                rhs_3 = velocity_3
                rhs_4 = G*M2*(x2 - x1)/|x2 - x1|^3 + G*M3*(x3 - x1)/|x3 - x1|^3
                rhs_5 = G*M3*(x3 - x2)/|x3 - x2|^3 + G*M1*(x1 - x2)/|x1 - x2|^3
                rhs_6 = G*M2*(x2 - x3)/|x2 - x3|^3 + G*M1*(x1 - x3)/|x1 - x3|^3
    initial parameters/time setup/initial conditions
    Calculate numerical solutions
    Comparison of numerical methods
"""

import matplotlib.animation as ani
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as spi

plt.close('all')

# differential equations set up
def rhs(time, state):
    
    # 18 equations: 3 2nd order * 2 for separating, * 3 for vector components
    position_1x = state[0]
    position_2x = state[1]
    position_3x = state[2]
    velocity_1x = state[3]
    velocity_2x = state[4]
    velocity_3x = state[5]
    position_1y = state[6]
    position_2y = state[7]
    position_3y = state[8]
    velocity_1y = state[9]
    velocity_2y = state[10]
    velocity_3y = state[11]
    position_1z = state[12]
    position_2z = state[13]
    position_3z = state[14]
    velocity_1z = state[15]
    velocity_2z = state[16]
    velocity_3z = state[17]
    
    position_1 = np.array([position_1x, position_1y, position_1z])
    position_2 = np.array([position_2x, position_2y, position_2z])
    position_3 = np.array([position_3x, position_3y, position_3z])
    
    # magnitudes of distances
    r1_2 = np.sqrt(np.sum((position_2 - position_1)**2))
    r1_3 = np.sqrt(np.sum((position_3 - position_1)**2))
    r2_3 = np.sqrt(np.sum((position_3 - position_2)**2))
    
    # right hand sides
    rhs = np.zeros(18)
    rhs[0] = velocity_1x
    rhs[1] = velocity_2x
    rhs[2] = velocity_3x
    rhs[3] = (G*M2*(position_2x - position_1x)/r1_2**3 + G*M3*(position_3x -
       position_1x)/r1_3**3)
    rhs[4] = (G*M3*(position_3x - position_2x)/r2_3**3 + G*M1*(position_1x -
       position_2x)/r1_2**3)
    rhs[5] = (G*M2*(position_2x - position_3x)/r2_3**3 + G*M1*(position_1x -
       position_3x)/r1_3**3)
    rhs[6] = velocity_1y
    rhs[7] = velocity_2y
    rhs[8] = velocity_3y
    rhs[9] = (G*M2*(position_2y - position_1y)/r1_2**3 + G*M3*(position_3y -
       position_1y)/r1_3**3)
    rhs[10] = (G*M3*(position_3y - position_2y)/r2_3**3 + G*M1*(position_1y -
       position_2y)/r1_2**3)
    rhs[11] = (G*M2*(position_2y - position_3y)/r2_3**3 + G*M1*(position_1y -
       position_3y)/r1_3**3)
    rhs[12] = velocity_1z
    rhs[13] = velocity_2z
    rhs[14] = velocity_3z
    rhs[15] = (G*M2*(position_2z - position_1z)/r1_2**3 + G*M3*(position_3z -
       position_1z)/r1_3**3)
    rhs[16] = (G*M3*(position_3z - position_2z)/r2_3**3 + G*M1*(position_1z -
       position_2z)/r1_2**3)
    rhs[17] = (G*M2*(position_2z - position_3z)/r2_3**3 + G*M1*(position_1z -
       position_3z)/r1_3**3)
    return rhs

def total_energy(state):
    # define positions between objects for potential energies
    position_1 = np.array([state[0], state[6], state[12]])
    position_2 = np.array([state[1], state[7], state[13]])
    position_3 = np.array([state[2], state[8], state[14]])
    r1_2 = np.sqrt(np.sum((position_2 - position_1)**2))
    r1_3 = np.sqrt(np.sum((position_3 - position_1)**2))
    r2_3 = np.sqrt(np.sum((position_3 - position_2)**2))
    kinetic = 0.5*(M1*(state[3]**2 + state[9]**2 + state[15]**2) +
               M2*(state[4]**2 + state[10]**2 + state[16]**2) +
               M3*(state[5]**2 + state[11]**2 + state[17]**2))
    potential = -G*(M1*M2/r1_2 + M1*M3/r1_3 + M2*M3/r2_3)
    return kinetic + potential

def total_momentum(state):
    v1 = np.array([state[3], state[9], state[15]])
    v2 = np.array([state[4], state[10], state[16]])
    v3 = np.array([state[5], state[11], state[17]])
    return M1*v1 + M2*v2 + M3*v3

G = 4*np.pi**2 # AU^2/(year^2 * Msun)
M1 = 1 # Solar masses
M2 = 1
M3 = 1

t_min = 0 # years
t_max = 500 # years
N = 10000 # fixed 100 steps so animation isn't too slow
times = np.linspace(t_min, t_max, N)

# Initial conditions; this set happens to remain bounded for 250+ years
# depending on the method used
init_cond = np.zeros(18)
init_cond[0] = 1 # r1x; Astronomical units
init_cond[1] = -0.5 # r2x
init_cond[2] = -0.5 # r3x
init_cond[3] = 0 + 1 # v1x
init_cond[4] = -3.3*np.sqrt(3) + 1 # v2x
init_cond[5] = 3.3*np.sqrt(3) + 1 # v3x
init_cond[7] = np.sqrt(3)/2 # r2y
init_cond[8] = -np.sqrt(3)/2 # r3y
init_cond[9] = 6.6 # v1y
init_cond[10] = -3.3 # v2y
init_cond[11] = -3.3 # v3y
init_cond[15] = 1 # v1z (setting all to 1 changes the results interestingly)
init_cond[16] = 1 # v2z
init_cond[17] = 1 # v3z

'''
Method options:
    RK45
    RK23
    Radau
    BDF
    LSODA
'''

# Methods chosen were the most interesting and remained bounded for a 500
# year interval. RK23 notably loses and then regains a lot of energy around 230
# years.

method1 = 'RK23'
method2 = 'LSODA'

solution1 = spi.solve_ivp(rhs, [t_min, t_max], init_cond, t_eval = times,
                          method = method1)
solution2 = spi.solve_ivp(rhs, [t_min, t_max], init_cond, t_eval = times,
                          method = method2)

# define energy and momentum components for both methods
energy = np.zeros((len(times), 2))
momentum = np.zeros((len(times), 2, 3))
for i in range(N):
    energy[i, 0] = total_energy(solution1['y'][:, i])
    momentum[i, 0] = total_momentum(solution1['y'][:, i])
    energy[i, 1] = total_energy(solution2['y'][:, i])
    momentum[i, 1] = total_momentum(solution2['y'][:, i])

fig = plt.figure()
# initialize points to represent bodies and text to be animated
line1, = plt.plot([], [], 'ok')
line2, = plt.plot([], [], 'or')
time_text = plt.text(0.2, 0.2, '', transform=plt.gcf().transFigure)

# initialization function for animation
def init():
    line1.set_data([], [])
    line2.set_data([], [])
    time_text.set_text('')
    return line1, line2, time_text
# what to animate at every frame
def animate(i):
    # scaling speeds up the animation
    scaling = 10
    i*=scaling
    # plot the trajectories at given frame
    x1 = [solution1['y'][0][i], solution1['y'][1][i], solution1['y'][2][i]]
    y1 = [solution1['y'][6][i], solution1['y'][7][i], solution1['y'][8][i]]
    x2 = [solution2['y'][0][i], solution2['y'][1][i], solution2['y'][2][i]]
    y2 = [solution2['y'][6][i], solution2['y'][7][i], solution2['y'][8][i]]
    line1.set_data(x1,y1)
    line2.set_data(x2,y2)
    # change the text to reflect current age
    time_text.set_text('years: ' + str(int(i*t_max/(1000*scaling))))
    return line2, time_text#line2, time_text

# plot trajectories
# One note is that going past 100 years the bodies will eventually
# diverge for all methods, but how long it will take for them to do so and
# their resulting trajectories varies significantly. They are described below
# for the various methods. Results are using only 1000 points for various
# timescales
# RK45: ~250 years: All 3 orbits are unbound
# RK23: ~250 years: 2 bodies "crash" into each other but orbit each other
#       rapidly rather than diverge, contnue orbits with 3rd body. Orbits
#       become unbound at about 550 years, where 2 bodies are still rapidly
#       orbiting each other.
# Radau: ~150 years: One orbit becomes unbound while the other 2 orbit each
#       other stablely
# BDF: ~120 years: Similar to Radau but the direction of the unbound orbits
#       change, and also one of the bodies is ejected at a much higher velocity
# LSODA: ~250 years: 2 orbits are no longer stable, and what happens after
#       is timestep dependent. 500 year timescale best represents a completely
#       nonperiodic, yet bounded, orbit
#plt.plot(solution1['y'][0], solution1['y'][6], '--k', label =
#         'trajectories using method = ' + method1)
#plt.plot(solution1['y'][1], solution1['y'][7], '--k')
#plt.plot(solution1['y'][2], solution1['y'][8], '--k')
plt.plot(solution2['y'][0], solution2['y'][6], '--r', label =
         'trajectories using method = ' + method2)
plt.plot(solution2['y'][1], solution2['y'][7], '--r')
plt.plot(solution2['y'][2], solution2['y'][8], '--r')
plt.title("Trajectories")
plt.xlabel("x distance (AU)")
plt.ylabel("y distance (AU)")
plt.legend()
# implement animation
ani.FuncAnimation(fig, animate, init_func=init, frames=1000, interval=10,
                  blit = True)

# Energy is conserved adequately while the bodies are far away from each
# other but while they are close, the energy level of the system spikes
# up or down for all methods. The methods are listed from best to worst in
# terms of how well they conserve energy for the given trajectories for 100
# year trials:
# Radau: Conserves energy very nicely
# RK45 & RK23: loses energy
# LSODA: loses slightly more energy than RK23
# BDF: loses substantially more energy
# An interesting note is that RK45 both roughly lose energy at the same rate
# but alternate how well they are doing over time. The trajectories also split
# apart and the reconverge
plt.figure()
plt.plot(times, energy[:, 0], label = "method = " + method1)
plt.plot(times, energy[:, 1], label = "method = " + method2)
plt.title("Energy over time")
plt.xlabel("time (AU)")
plt.ylabel("Energy (Msun*AU^2*year^-2)")
plt.legend()

# momentum seems to be conserved relatively well for all methods
# X momentum (100 year trials)
# RK45, RK23, LSODA: conserves perfectly
# Radau: gains momentum (order e-12)
# BDF: momentum varies, order e-11
plt.figure()
plt.plot(times, momentum[:, 0, 0], label = 'method = ' + method1)
plt.plot(times, momentum[:, 1, 0], label = 'method = ' + method2)
plt.title("X momentum over time")
plt.xlabel("time (AU)")
plt.ylabel("momentum (Msun*AU*year^-1)")
plt.legend()

# Y momentum (100 year trials)
# RK23: gains momentum in +y, order e-15
# RK45: gains momentum in -y, order e-15
# LSODA: gains momentum in +y, order e-14
# Radau: gains momentum in +y, order e-9
# BDF: loses momentum in +y, order e-9
plt.figure()
plt.plot(times, momentum[:, 0, 1], label = 'method = ' + method1)
plt.plot(times, momentum[:, 1, 1], label = 'method = ' + method2)
plt.title("Y momentum over time")
plt.xlabel("time (AU)")
plt.ylabel("momentum (Msun*AU*year^-1)")
plt.legend()

# Z momentum
# This is conserved perfectly for all methods, as is expected
plt.figure()
plt.plot(times, momentum[:, 0, 2], label = 'method = ' + method1)
plt.plot(times, momentum[:, 1, 2], label = 'method = ' + method2)
plt.title("Z momentum over time")
plt.xlabel("time (AU)")
plt.ylabel("momentum (Msun*AU*year^-1)")
plt.legend()
