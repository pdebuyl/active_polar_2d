#!/usr/bin/env python3
"""Execute a Langevin simulation"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int)
parser.add_argument('--steps', type=int, default=2**11)
parser.add_argument('--gamma', type=float, default=1)
parser.add_argument('-T', type=float, default=1)
parser.add_argument('--v0', type=float, default=0)
parser.add_argument('--gamma-par', type=float, default=1)
parser.add_argument('--gamma-perp', type=float, default=1)
parser.add_argument('--gamma-r', type=float, default=1)
parser.add_argument('--gravity', type=float, default=0)
parser.add_argument('-L', type=float, default=10)
args = parser.parse_args()

import numpy as np
import matplotlib.pyplot as plt
import active_polar_2d
import scipy.signal
import math

if args.seed is not None:
    seed = args.seed
else:
    import time
    seed = int((time.time()*1000) % 2**32)

# Run the simulation

interval = 100
dt = 0.001

r_cut = 3**(1/6)
L = args.L


def wall_f(x):
    if x < r_cut:
        r = x
        sign = 1
    elif x > L-r_cut:
        r = L-x
        sign = -1
    else:
        return 0
    return sign*3*(3*r**-10 - r**-4)


grav = float(args.gravity)
bond_l = 2


def f(x, y, theta):
    c_th = math.cos(theta)
    s_th = math.sin(theta)
    x1 = x + c_th
    f1 = wall_f(x1)
    x2 = x - c_th
    f2 = wall_f(x2)
    return f1 + f2 - grav, 0, s_th*bond_l/2*(f1-f2)


x, y, theta = \
active_polar_2d.integrate_OD_2d_theta(L/2, 0, 0,
                                      args.gamma_par, args.gamma_perp,
                                      args.gamma_r, args.T, args.v0, dt,
                                      interval, args.steps, f=f, seed=seed)

# Display the results
N = len(x)
t = dt*interval*np.arange(N)

plt.figure()
ax1 = plt.subplot(311)
plt.plot(y, x)
plt.ylabel(r'$x$')

ax2 = plt.subplot(312)
plt.plot(t, np.mod(theta, 2*np.pi)/np.pi)
plt.ylabel(r'$\theta/\pi$')

ax3 = plt.subplot(325)
plt.hist(x, normed=True, bins=32)

ax4 = plt.subplot(326)
plt.hist(np.mod(theta, 2*np.pi), normed=True, bins=32)

plt.show()
