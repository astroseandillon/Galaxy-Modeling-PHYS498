# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 23:20:27 2020

@author: Seam Dillon

This file is to create an object that we can refer to our gravitational bodies.
I figure that if we can describe each of galaxies by their conditions, which 
can help simplify our code.
"""

class Galaxy():
    """
    This is just some ideas for properties we can give to our galaxies. Anyone 
    who has any idea of other things to add, feel free to add them and we can 
    figure it out together.
    
    As it stands, we begin by inserting in all of our properties of our galaxy
    into the initial state 
    
    """
    
    def __init__(self, stellar_mass, radius, phi_min, phi_max, theta_min, theta_max, barred=False):
        self.stellar_mass = stellar_mass
        self.radius = radius
        self.phi_min = phi_min
        self.phi_max = phi_max
        self.theta_min = theta_min
        self.theta_max = theta_max
        return 

    def state(position, velocity):
        self.position[0] = x
        










