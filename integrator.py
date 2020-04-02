# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 10:13:25 2020

@author: Sean Dillon, Bjorn Larsen

Resources:
https://pdfs.semanticscholar.org/2103/85ca3ad4c53f946c7dc0a0b259cd38b6bdf7.pdf
https://en.wikipedia.org/wiki/Symplectic_integrator
Candy & Rozmus: A Symplectic Integration Algorithm for Separable Hamiltonian Functions

q = position coordinates
p = momentum coordinates
H = hamiltonian
z = (q,p) = canonical coordinates

H(p,q) = T(p) + V(q)

T = kinetic energy
V = potential energy

Hamilton's equations:
F = dp/dt = -∂H/∂q
P = dq/dt = ∂H/∂p

The integration scheme:
    Do for i in range(4):
        p[i] = p[i] + b[i]F(q[i-1])*dt
        q[i] = q[i-1] + a[i]P(p[i])
    Integrated variables: (q[3], p[3]) at t = t0 + dt
"""

import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

# coefficients for fourth-order integration
a = np.array([
    (2 + 2**(1/3) + 2**(-1/3))/6,   # a1
    (1 - 2**(1/3) - 2**(-1/3))/6,   # a2
    (1 - 2**(1/3) - 2**(-1/3))/6,   # a3
    (2 + 2**(1/3) + 2**(-1/3))/6    # a4
])

b = np.array([
    0,                      # b1
    (2 - 2**(1/3))**-1,     # b2
    (1 - 2**(2/3))**-1,     # b3
    (2 - 2**(1/3))**-1      # b4
])

# time parameters
t_min = 0
t_max = 10
Nt = 10000
times = np.linspace(t_min, t_max, Nt)
dt = times[1]-times[0]

# initial conditions
q0 = np.pi/2
p0 = 0
M = 1
g = 9.8
L = 1


def V(q):
    """
    POTENTIAL ENERGY
    """
    return (g/L)*(1-np.cos(q))

def F(q):
    return -(g/L)*np.sin(q)

def T(p):
    """
    KINETIC ENERGY
    """
    return (p**2)/(2*M)

def P(p):
    return p/M


def H(coordinates):
        """The Hamiltonian is the sum of kinetic and potential energy."""
        q = coordinates[0,:]
        p = coordinates[1,:]
        return T(p) + V(q)

def first_order_equation(t, q0, p0):
    q = np.zeros(Nt)
    p = np.zeros(Nt)
    p[0] = p0
    q[0] = q0
    for i in range(1,Nt):
        p[i] = p[i-1] + dt*F(q[i-1])
        q[i] = q[i-1] + dt*P(p[i])
    z = np.array((q, p))
    return z

def fourth_order_equation(t, q0, p0):
    q = np.zeros(Nt)
    p = np.zeros(Nt)
    p[0] = p0
    q[0] = q0
    for i in range(1,Nt):
        p_temp = p[i-1]
        q_temp = q[i-1]
        for j in range(4):
            p_temp += b[j]*F(q_temp)*dt
            q_temp += a[j]*P(p_temp)*dt
        p[i] = p_temp
        q[i] = q_temp
    z = np.array((q, p))
    return z

testfun = first_order_equation(times, q0, p0)
new_integrator = fourth_order_equation(times, q0, p0)

x = testfun[0,:]
y = testfun[1,:]

new_x = new_integrator[0,:]
new_y = new_integrator[1,:]

plt.figure()
plt.plot(times, x, 'o', label="first order")
plt.plot(times, new_x, '.', label="fourth order")
plt.title("Symplectic Integrator Test")
plt.ylabel("Position")
plt.xlabel("Time")
plt.legend()

plt.figure()
plt.plot(times, H(testfun), label="first order")
plt.plot(times, H(new_integrator), label="fourth order")
plt.title("conservation of energy i think")
plt.xlabel('time')
plt.ylabel("hamiltonian")
plt.legend()


plt.show()





