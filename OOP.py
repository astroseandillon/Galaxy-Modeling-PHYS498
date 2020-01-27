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
    
    def speak(self, sound):
        return "{} says {}".format(self.name, sound)
        
    def description(self):
        return "{} is {} years old".format(self.name, self.age)    
        
        
        
        
class Mutt(Dog):
    
    def run(self, speed):
        return "{} runs {}".format(self.name, speed)     
        
 



      

molly = Mutt('Molly', 10)
print(molly.description()) 

print(molly.run('slowly'))
 
