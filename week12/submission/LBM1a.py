"""
Applied class (starter): flow past a cylinder.

Velocity inlet (Zou/He, west), extrapolated outlet (east), slip walls
(north/south) and a bounce-back cylinder, D = 10 nodes, Re = 100.
Look for the TODO comments for where to add the pressure probe.
"""
from time import time
from git import Repo
import numpy as np
import matplotlib.pyplot as plt


# Lattice constants (D2Q9)
DISCRETE_VELOCITIES = 9
DIRECTIONS = np.array([(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1),
                       (1, 1), (-1, 1), (-1, -1), (1, -1)]).T
BOUNCED = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])
FULLSLIP = np.array([0, 1, 4, 3, 2, 8, 7, 6, 5])  # reflect the y-component only
WEIGHTS = np.array([4/9,
                    1/9, 1/9, 1/9, 1/9,
                    1/36, 1/36, 1/36, 1/36], dtype=np.float64)
CS2 = 1/3  # lattice speed of sound squared


def f_equilibrium(rho, u_x, u_y):
    feq = np.zeros((DISCRETE_VELOCITIES,) + rho.shape)
    for i in range(DISCRETE_VELOCITIES):
        cu = DIRECTIONS[0, i]*u_x + DIRECTIONS[1, i]*u_y
        feq[i] = WEIGHTS[i]*rho*(1 + cu/CS2 + cu**2/(2*CS2**2) - (u_x**2 + u_y**2)/(2*CS2))
    return feq


def moments(f):
    rho = np.sum(f, axis=0)
    u_x = np.einsum("ijk, i -> jk", f, DIRECTIONS[0, :]) / rho
    u_y = np.einsum("ijk, i -> jk", f, DIRECTIONS[1, :]) / rho
    return rho, u_x, u_y


def collisionMRT(f, S2, S3=0, S4=0):
    # Moment-space collision: each moment R relaxes as R -> S*(R - R_eq) + R_eq.
    # S2 = 1 - omega for the stress moments (sets the viscosity); S3 = S4 = 0
    # puts the higher (non-hydrodynamic) moments straight to equilibrium.
    f_post = np.zeros_like(f)
    rho = f[8] + f[7] + f[6] + f[5] + f[4] + f[3] + f[2] + f[1] + f[0]
    Jx = f[8] - f[7] - f[6] + f[5] - f[3] + f[1]
    Jy = -f[8] - f[7] + f[6] + f[5] - f[4] + f[2]
    R3 = -f[4] - f[2] - f[0] + ( f[8] + f[7] + f[6] + f[5] + f[3] + f[1] )*2.
    R4 = -f[3] - f[1] - f[0] + ( f[8] + f[7] + f[6] + f[5] + f[4] + f[2] )*2.
    R5 = ( -f[8] + f[7] - f[6] + f[5] )*3.
    R6 = f[4] - f[2] + ( -f[8] - f[7] + f[6] + f[5] )*2.
    R7 = f[3] - f[1] + ( f[8] - f[7] - f[6] + f[5] )*2.
    R8 = f[0] + ( -f[4] - f[3] - f[2] - f[1] + ( f[8] + f[7] + f[6] + f[5] )*2. )*2.

    R3 = S2*(R3 - Jx*Jx/rho*3.) + Jx*Jx/rho*3.
    R4 = S2*(R4 - Jy*Jy/rho*3.) + Jy*Jy/rho*3.
    R5 = S2*(R5 - Jy*Jx/rho*3.) + Jy*Jx/rho*3.
    R6 = S3*R6
    R7 = S3*R7
    R8 = S4*R8

    f_post[0] = ( R8 + ( -R4 - R3 + rho*2. )*2. )/9.
    f_post[1] = ( -R8 - R4 - R7*3. + ( R3 + rho + Jx*3. )*2. )/18.
    f_post[2] = ( -R8 - R3 - R6*3. + ( R4 + rho + Jy*3. )*2. )/18.
    f_post[3] = ( -R8 - R4 + R7*3. + ( R3 + rho - Jx*3. )*2. )/18.
    f_post[4] = ( -R8 - R3 + R6*3. + ( R4 + rho - Jy*3. )*2. )/18.
    f_post[5] = ( R8 + R4 + R3 + rho + ( R7 + R6 + R5 + Jy + Jx )*3. )/36.
    f_post[6] = ( R8 + R4 + R3 + rho + ( -R7 + R6 - R5 + Jy - Jx )*3. )/36.
    f_post[7] = ( R8 + R4 + R3 + rho + ( -R7 - R6 + R5 - Jy - Jx )*3. )/36.
    f_post[8] = ( R8 + R4 + R3 + rho + ( R7 - R6 - R5 - Jy + Jx )*3. )/36.
    return f_post


cyl_r = 5
NX = cyl_r*2 * 25; NY = cyl_r*2 * 20  # L x H = 25D x 20D
cyl_cx = NX // 4; cyl_cy = NY // 2

N_ITER = 20_000
PLOT_EVERY = N_ITER // 100

RE = 100
U0 = 0.05
collision = 'MRT'  # 'BGK'

nu    = U0 * 2 * cyl_r / RE
tau   = 3 * nu + 0.5
omega = 1/tau
print(f"tau = {tau:.4f}, omega = {omega:.4f}, Ma = {U0/np.sqrt(CS2):.3f}")

# TODO: define the pressure probe location and a list to store its history

# Solid flags: the cylinder
xx, yy = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
solid = np.zeros((NX, NY), dtype=bool)
solid[(xx - cyl_cx)**2 + (yy - cyl_cy)**2 <= cyl_r**2] = True

# Initialise with uniform flow U0, slightly faster in the lower half to trigger shedding
rho = np.ones((NX, NY))
u_x = U0*np.ones((NX, NY)); u_y = np.zeros((NX, NY))
u_x[:, 2:NY//2] = 1.1*U0
u_x[solid] = 0
f_prec = f_equilibrium(rho, u_x, u_y)

fig, ax = plt.subplots(figsize=(10, 10*NY/NX))  # TODO: add a second subplot for the probe pressure
tic = time()
for iter in range(N_ITER):
    #1) Macroscopic quantities
    rho, u_x, u_y = moments(f_prec)
    # TODO: record the probe pressure, p = rho c_s^2

    if iter % PLOT_EVERY == 0:
        toc = time()
        print(f"{iter=}: fluid mass = {np.sum(rho[~solid]):.3f}, max |u| = {np.max(np.sqrt(u_x**2 + u_y**2)):.4f}"
              f" ... {NX*NY*PLOT_EVERY/(toc - tic)/1e6:.2f} MLUPS")
        ax.clear()
        ax.imshow(np.sqrt(u_x**2 + u_y**2).T, origin='lower', cmap='viridis')
        ax.add_patch(plt.Circle((cyl_cx, cyl_cy), cyl_r, color='red'))
        ax.set_title(f'Iteration {iter}')
        plt.pause(0.001)
        tic = time()

    #2) Equilibrium
    f_eq = f_equilibrium(rho, u_x, u_y)

    #3) Collision; solid nodes do not collide
    f_solid = f_prec[:, solid]
    if collision == 'BGK':
        f_post = f_prec - omega*(f_prec - f_eq)
    elif collision == 'MRT':
        f_post = collisionMRT(f_prec, S2=1 - omega)
    f_post[:, solid] = f_solid

    #4) Streaming
    for i in range(DISCRETE_VELOCITIES):
        f_prec[i, :, :] = np.roll(f_post[i, :, :], DIRECTIONS[:, i], axis=(0, 1))

    #5) Bounce-back: reverse every population that arrived at a solid node
    f_tmp = f_prec.copy()
    for i in range(DISCRETE_VELOCITIES):
        f_prec[i, solid] = f_tmp[BOUNCED[i], solid]

    #5a) Slip walls (north/south): reflect the y-component, keep the x-component
    for i in (4, 7, 8):
        f_prec[i, :, -1] = f_prec[FULLSLIP[i], :, -1]
    for i in (2, 5, 6):
        f_prec[i, :, 0] = f_prec[FULLSLIP[i], :, 0]

    #5b) Inlet (west): Zou/He velocity boundary, unknowns f1, f5, f8
    rho_w = 1/(1-U0) * (f_prec[0,0] + f_prec[2,0] + f_prec[4,0] + 2*(f_prec[3,0] + f_prec[6,0] + f_prec[7,0]))
    f_prec[1,0] =  f_prec[3,0] + 2/3*rho_w*U0
    f_prec[5,0] = -0.5*f_prec[2,0] + 0.5*f_prec[4,0] + f_prec[7,0] + 1/6*rho_w*U0
    f_prec[8,0] =  0.5*f_prec[2,0] - 0.5*f_prec[4,0] + f_prec[6,0] + 1/6*rho_w*U0

    #5c) Outlet (east): copy rho and u from the column inside, set f = f_eq
    #    (simpler alternative: zero-gradient for f3, f6, f7, i.e. f_prec[i,-1] = f_prec[i,-2])
    rho_out, u_x_out, u_y_out = moments(f_prec[:, -2:-1, :])
    f_prec[:, -1:, :] = f_equilibrium(rho_out, u_x_out, u_y_out)

# TODO: save the pressure trace to a csv file and save the figure
plt.show()
