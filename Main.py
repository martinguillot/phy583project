# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 15:00:52 2019

@author: marti
"""

from Orbit import *
from Binary import *
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation

N = 2000000
orbit1 = Orbit.Orbit(N, 0.01, 1, [20, 0.], [0.1, 0.1])
orbit2 = Orbit.Orbit(N, 0.1, 1., [-20, 0.], [-0.1, -0.1])
binarytest = Binary(orbit1, orbit2)

binarytest.compute()
binarytest.visualize()
plt.show()

fig = plt.figure()
ax = plt.axes(xlim=(0, 2), ylim=(-2, 2))
line, = ax.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    x = binarytest.orbi1.grid[:i][0][0]
    y = binarytest.orbi1.grid[:i][0][1]
    line.set_data(x, y)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=2000, interval=20, blit=True)

plt.show()
