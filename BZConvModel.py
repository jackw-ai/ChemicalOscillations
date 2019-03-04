'''
Simulated equations:
    a_{t+1} = a_t + a_t(alpha b_t - gamma c_t)
    
    b_{t+1} = b_t + b_t(beta c_t - alpha a_t)
    
    c_{t+1} = c_t + c_t(gamma a_t - beta  b_t)
'''

import numpy as np
from scipy.signal import convolve2d
import matplotlib.pyplot as plt
from matplotlib import animation

from IPython.display import HTML

WIDTH, HEIGHT = 600, 400

class BZSimConv():
    ''' Improved BZ model that adds parameter constants,
        uses numpy for improved efficiency,
        and uses scipy convolutions to find average
    '''
    
    def __init__(self, width, height, alpha = 1.2, beta = 1, gamma = 1):
        # width, height of the image
        self.w, self.h = width, height
    
        # reaction constants
        self.alpha, self.beta, self.gamma = alpha, beta, gamma
    
        # initialize grid with random amounts of a, b and c
        self.data = np.random.random(size=(2, 3, self.h, self.w))
        
        # bit mask
        self.p = 0
        
    def reset(self, alpha = 1.2, beta = 1, gamma = 1):
        ''' resets the grid matrices and updates reaction constants '''
        
        # re-initialize grid with random amounts of a, b and c
        self.data = np.random.random(size=(2, 3, self.h, self.w))

        # update reaction constants
        self.alpha, self.beta, self.gamma = alpha, beta, gamma
        
    def step(self, p = -1):
        """ updates data[p] using data[q] for a time step """
    
        # we pass in a separate p for easier animating
        p = p if p > -1 else self.p
    
        # get next grid
        q = 1 - p
    
        s = np.zeros((3, self.h, self.w))
    
        m = np.ones((3, 3)) / 9
    
        # count the average amount of each species in the 9 cells around each cell
        for k in range(3):
            # uses convolution with 3x3 matrix m for faster computation
            s[k] = convolve2d(self.data[p,k], m, mode = 'same', boundary = 'wrap')

        # reaction equations
        self.data[q, 0] = s[0] + s[0]*(self.alpha * s[1] - self.gamma * s[2])
        self.data[q, 1] = s[1] + s[1]*(self.beta * s[2] - self.alpha * s[0])
        self.data[q, 2] = s[2] + s[2]*(self.gamma * s[0] - self.beta * s[1])

        # keeps species concentrations within [0, 1]
        np.clip(self.data[q], 0, 1, self.data[q])
    
    
    def simulate(self, save = False, colormap = plt.cm.winter):
        ''' generates video of simulation '''
        
        # Set up the image
        fig, ax = plt.subplots()
        im = ax.imshow(self.data[0,0], cmap = colormap)
        ax.axis('off')

        def animate(i):
            """ generate frame i for plt animation """
            arr = self.step(i % 2)
            im.set_array(self.data[i % 2, 0])
            return [im]
        
        anim = animation.FuncAnimation(fig, animate, frames = 200, interval = 100, blit=False)

        if save:
            anim.save(filename = 'bz.mp4', fps = 15)
            
        return HTML(anim.to_html5_video())


conv = BZSimConv(WIDTH, HEIGHT)
conv.simulate()

conv.reset(1, 1, 1)
conv.simulate(save = True, colormap = 'jet')

conv.reset(1.2, 1, 1)
conv.simulate(save = True)

conv.reset(1, 1, 1.4)
conv.simulate()

conv.reset(1, 1.3, 1)
conv.simulate(save = True)
