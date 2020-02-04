# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 23:20:27 2020

@author: Seam Dillon

This file is to create an object that we can refer to our gravitational bodies.
I figure that if we can describe each of galaxies by their conditions, which 
can help simplify our code.
"""
import numpy as np
class Galaxy():
    """
    This is just some ideas for properties we can give to our galaxies. Anyone 
    who has any idea of other things to add, feel free to add them and we can 
    figure it out together.
    
    CREATE CHILD CLASSES
    
    """
    
    def __init__(self, stellar_mass, radius, phi_min, phi_max, theta_min, theta_max):
        self.stellar_mass = stellar_mass
        self.radius = radius
        self.phi_min = phi_min
        self.phi_max = phi_max
        self.theta_min = theta_min
        self.theta_max = theta_max
        return 

m51 = Galaxy(10, 9, 8, 7, 6, 5)

print(m51.stellar_mass)
print(m51.radius)
print(m51.phi_min)
print(m51.phi_max)
print(m51.theta_min)
print(m51.theta_max)
