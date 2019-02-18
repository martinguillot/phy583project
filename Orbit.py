# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 00:26:44 2019

@author: marti
"""
import numpy as np
import matplotlib.pyplot as plt
#G= 6.67408*10^(-11)
G = 1

class Orbit:
    "contain data about a massive body orbit"
    def __init__(self, N, dt, M, Pinit, Vinit):
        "Create an empty orbit object with an initial speed and position V0 and P0"
        self.N = N
        self.dt = dt
        self.M = M
        self.grid = np.zeros((N,2,2))
        self.grid[0][0] = Pinit
        self.grid[0][1] = Vinit
    
    def energy(self, orbit2, t):
        "gives the energy of self object at a time increment t"
        r = ((self.grid[t][0][0]-orbit2[t][0][0])^2+(self.grid[t][0][1]-orbit2[t][0][1])^2)^(1/2)
        return 0.5*self.M*(self.grid[t][1][0]^2+self.grid[t][1][1]^2)-self.M*G*orbit2.M/r
    
    def kinetic(self, t):
        "return the kinetic enery of the system at a time t"
        return 0.5*self.M*(self.grid[t][1][0]^2+self.grid[t][1][1]^2)
    
    def acceleration(self):
        A = np.zeros((2,self.N))
        A = (self.grid[1:][1]-self.grid[:-2][1])/self.dt
        A[-1] = A[-2]
        return A
    
    def visualize(self):
        X = np.zeros(self.N)
        Y = np.zeros(self.N)
        for i in range(self.N):
            X[i] = self.grid[i][0][0]
            Y[i] = self.grid[i][0][1]
        plt.plot(X, Y)