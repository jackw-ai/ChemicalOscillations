# Generates a visual representation of a Belousov - Zhabotinsky Reaction

''' reaction equations
src : http://discovery.ucl.ac.uk/17241/1/17241.pdf

Simulated equations:
a_t+1 = a_t + a_t(b_t - c_t)
b_t+1 = b_t + b_t(c_t - a_t)
c_t+1 = c_t + c_t(a_t - b_t)
'''

import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

# window dimensions
WIDTH, HEIGHT = 400, 400

class BZSim():
    def __init__(self, width, height):
        self.w = width
        self.h = height

        # initialize grid representation
        ## since we need the previous instance, we store two matrices
        ## and switch between then at each timestep
        self.a = np.zeros((2, WIDTH, HEIGHT))
        self.b = np.zeros((2, WIDTH, HEIGHT))
        self.c = np.zeros((2, WIDTH, HEIGHT))

        # bit index representing which matrix is curr timestep
        self.p = 0
        
    def setup(self):
        ''' initializes a, b, c matrices where elements ~ unif(0,1) '''
        for i in range(self.w):
            for j in range(self.h):
                self.a[self.p][i][j] = random.random()
                self.b[self.p][i][j] = random.random()
                self.c[self.p][i][j] = random.random()

    @staticmethod
    def clamp(n, small, big):
        '''  helper, enforces that small <= n <= big '''
        return max(min(big, n), small)
    
    def step(self):
        ''' one step in the simulation '''
        
        # inverse of the timestep bit
        q = 1 - self.p
        
        for i in range(self.w):
            for j in range(self.h):
                c_a, c_b, c_c = self.get_adjacents(i, j)

                # reaction equations
                self.a[q][i][j] = self.clamp(c_a + c_a * (c_b - c_c), 0, 1)
                self.b[q][i][j] = self.clamp(c_b + c_b * (c_c - c_a), 0, 1)
                self.c[q][i][j] = self.clamp(c_c + c_c * (c_a - c_b), 0, 1)

        # inverts the bit since curr timestep is now prev timestep
        self.p = 1 - self.p

    def get_adjacents(self, i, j):
        ''' returns the average value of the surrounding entries '''
        c_a = c_b = c_c = 0.0

        for k in range(3):
            for l in range(3):

                # cheating alittle in edge case by just wrapping around matrix
                x = (i + k) % self.w
                y = (j + l) % self.h
                
                c_a += self.a[self.p][x][y]
                c_b += self.b[self.p][x][y]
                c_c += self.c[self.p][x][y]

        # average across the nine values
        return c_a / 9, c_b / 9, c_c / 9

    def animate(self, n_frames = 10, cmap = 'jet'):
        ''' animates grid using plt for n_frames using cmap as colormap'''
        
        fig = plt.figure()
        im = plt.imshow(self.a[self.p])

        ims = []
        ims.append([im])
        
        # add the next frames
        for _ in range(n_frames):
            self.step()
            im = plt.imshow(self.a[1 - self.p], cmap = cmap, animated = True)
            ims.append([im])

        anim = animation.ArtistAnimation(fig, ims, interval=100, blit=True, repeat_delay=2000)

        plt.show()

    def display(self, cmap = 'jet', initial = False):
        ''' displays current state '''
        p = self.p if initial else 1 - self.p
        
        plt.imshow(self.a[p], cmap = cmap)
        plt.show()
        
sim = BZSim(WIDTH, HEIGHT)

sim.setup()
sim.display(initial = True)
sim.step()
sim.display()
sim.step()
sim.display()

#sim.animate()
