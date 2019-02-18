# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 11:26:29 2019

@author: marti
"""

import Orbit
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

G = 1
#G= 6.67408*10^(-11)

class Binary:
    "Contains all the functions used to compute the evolution of a binary system"
    
    def __init__(self, orbit1, orbit2):
        "Takes two orbit-class object and creates a binary class object"
        self.orbit1 = orbit1
        self.orbit2 = orbit2
        self.M = orbit1.M+orbit2.M
        self.mu = orbit1.M*orbit2.M/self.M
        self.dt = orbit1.dt
        self.N = orbit1.N
    
    def distance(self, t):
        "Compute the distance between bodies at a time t"
        r = ((self.orbit1.grid[t][0][0]-self.orbit2.grid[t][0][0])**2+(self.orbit1.grid[t][0][1]-self.orbit2.grid[t][0][1])**2)**(0.5)
        return r
        
    def kineticenergy(self, t):
        "return the total kinetic energy of the system"
        return self.orbit1.kinetic(t)+self.orbit2.kinetic(t)
    
    def potentialenergy(self, t):
        return -2*G*self.mu*self.M/self.distance(t)
        
    def totalenergy(self, t):
        "compute the total energy of the system at a time t"
        return self.kineticenergy(t)+self.potentialenergy(t)
    
    def force(self, t):
        "return a 2x2 array of the force imposed on the system at a time t"
        force = np.zeros(2)
        u = [self.orbit1.grid[t][0][0]-self.orbit2.grid[t][0][0], self.orbit1.grid[t][0][1]-self.orbit2.grid[t][0][1]]
        force[0] = (-G*self.mu*self.M/self.distance(t)**3.)*u[0]
        force[1] = (-G*self.mu*self.M/self.distance(t)**3.)*u[1]
        return force
        
    def propagate(self, t):
        F = self.force(t)
        self.orbit1.grid[t+1][1] = self.orbit1.grid[t][1]+self.dt*F/self.orbit1.M
        self.orbit2.grid[t+1][1] = self.orbit2.grid[t][1]-self.dt*F/self.orbit2.M
        self.orbit1.grid[t+1][0] = self.orbit1.grid[t][0]+self.dt*self.orbit1.grid[t][1]
        self.orbit2.grid[t+1][0] = self.orbit2.grid[t][0]+self.dt*self.orbit2.grid[t][1]
    
    def compute(self):
        "takes an empty binary object with just the initial condition and compute the entire trajectories"
        for i in range(self.N-1):
            self.propagate(i)
        self.orbit1.grid[-1] = self.orbit1.grid[-2]
        self.orbit2.grid[-1] = self.orbit2.grid[-2]
    
    def visualize(self):
        self.orbit1.visualize()
        self.orbit2.visualize()
        
       
    