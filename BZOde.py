# ODE for B-Z Model
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# BZ constants
q, f, eps, delta = 3.1746e-5, 1, 0.0099, 2.4802e-5

# initial values
x0, y0, z0 = 0, 0, 0.1

# end time
t_e = 50 

# number of timesteps
n_iter = 10000

# t values
t = np.linspace(0, t_e, n_iter)

# reaction ODE
def bz_reaction(triple, t, q, f, eps, delta):
    x, y, z = triple # (x, y, z) tuple
    dx = (q*y - x*y + x*(1-x)) / eps
    dy = (-q*y - x*y + f*z) / delta
    dz = x-z
    return (dx, dy, dz)

f = odeint(bz_reaction, (x0, y0, z0), t, args=((q, f, eps, delta)))

# get results
x, y, z = f.T

# plot results
fig = plt.figure(figsize=(10, 15), facecolor='white')
fig.subplots_adjust(wspace=0.1, hspace=0.1)

ax1 = fig.add_subplot(3, 1, 1)
ax1.set_ylabel('concentration of bromous acid', fontsize=12)

ax2 = fig.add_subplot(3, 1, 2)
ax2.set_ylabel('concentration of bromide ions', fontsize=12)

ax3 = fig.add_subplot(3, 1, 3)
ax3.set_ylabel('concentration of cerium ions', fontsize=12)

ax1.set_title('Belousov - Zhabotinsky Reaction', fontsize = 20)

ax1.plot(t, x, 'b-')
ax2.plot(t, y, 'r-')
ax3.plot(t, z, 'g-')

plt.show()
