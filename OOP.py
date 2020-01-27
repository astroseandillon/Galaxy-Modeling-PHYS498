# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 21:58:48 2020

@author: seand
"""

#Object Oriented Programming
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import astropy as ap




class Dog:
    
    species = 'mammal'
    
    def __init__(self, name, age):
        self.name = name
        self.age = age
        
        
        
        

molly = Dog('Molly', 10)

print("{} is {}".format(molly.name, molly.age))

if molly.species == "mammal":
    print("{} is a {}!".format(molly.name, molly.species))
    
