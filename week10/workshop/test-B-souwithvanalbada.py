# test-A.py

from math import sin, pi
from numpy import zeros

def test_f(x):
    if 2<= x <= 3:
        return 1.0
    else:
        return 0.0

def van_albada_limiter(r):
    return (r**2 + r)/(r**2 + 1.0)

def update_in_time(f, a, L, ncells, dt, nsteps):
    dx = float(L)/float(ncells)
    xs = zeros(ncells)
    us = zeros(ncells)
    fs = zeros(ncells+1)

    # Initialise positions and starting values
    for i in range(ncells):
        xs[i] = dx*(i+0.5)
        us[i] = f(xs[i])

    # Now begin timestepping
    for n in range(nsteps):
        # 1. Compute fluxes.
        # 1a. At all INTERNAL interfaces
        for i in range(2, ncells):
            if us[i-1] - us[i-2] != 0.0:    
                # van albada limiter
                r = (us[i]-us[i-1])/(us[i-1]-us[i-2])
                phi = van_albada_limiter(r)
            else:
                phi = 1.0
            # second order upwind scheme
            fs[i] = us[i-1] + 0.5*phi*(us[i-1] - us[i-2])
        # 1b. At the BOUNDARY interfaces
        fs[0] = 0.0; fs[-1] = a*us[-1]

        # 2. Compute new u value
        for i in range(ncells):
            us[i] = us[i] - (dt/dx)*(fs[i+1] - fs[i])
        
    return xs, us

dt = 0.01
nsteps = 10
L = 5.0
a = 1.0
ncells = 100
dx = L/ncells

# generate analytical solution
xa = [ dx*(i+0.5) for i in range(ncells) ]
t_final = dt*nsteps
ua = [ test_f(x - a*t_final) for x in xa ]

# generate numerical solution
# Call update_in_time
xn, un = update_in_time(test_f, a, L, ncells, dt, nsteps)
# Save in xn, un

import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.plot(xa, ua, "-", xn, un, "o")
plt.show()


