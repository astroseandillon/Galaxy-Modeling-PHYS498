# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 10:13:25 2020

@author: Sean Dillon

After finding out that i have no fucking idea of what the symplectic integrator
code does, I'm building one from scratch. I am a dumb bitch and we shall now 
witness the power of street knowledge

https://pdfs.semanticscholar.org/2103/85ca3ad4c53f946c7dc0a0b259cd38b6bdf7.pdf

https://en.wikipedia.org/wiki/Symplectic_integrator

q = position coordinates
p = momentum coordinates
H = hamiltonian
z = (q,p) = canonical coordinates

H(p,q) = T(p) + V(q)

T = kinetic energy
V = potential energy

dp/dt = - partial H/partial q

dq/dt = partial H/partial p

dz/dt = {z, H(z)} = D_h z

z(tau) = exp(tau D_h)z(0)

D_h dot = {dot, H}






"""
import numpy as np


c14 = 1/(2*(2-(2**(1/3))))

c_2_4 = (1-(2**(1/3)))/(2*(2-(2**(1/3))))

d_1_4 = 1/(2-(2**(1/3)))

d_2_4 = (-2**(1/3))*d_1_4

d_4_4 = 0.0

first_order = np.array([[1.0],
                        [1.0]])

second_order = np.array([[0.0, 1.0],
                         [0.5, 0.5]])

third_order = np.array([[1.0, -(2/3), (2/3)],
                         [-(1/24), 0.75, (7/24)]])

fourth_order= np.array([[c14, c_2_4, c_2_4, c14],
                        [d_1_4, d_2_4, d_1_4, d_4_4]])


timestep = np.arange(0, 100, 0.1)



    



   
def first_order_equation(t, M):
    q = np.zeros(1000)
    p = np.zeros(1000)
    p[0] = np.pi/2
    q[0] = 10
    for i in range(len(t)-1):
        q[i+1] = q[i] + first_order[0]*p[i+1]*t[i]/M
        p[i+1] = p[i] + first_order[1]*M*t[i]*p[i]*p[i]/q[i]
    z = np.array((q, p))
    return z



testfun = first_order_equation(timestep, 1.0)















