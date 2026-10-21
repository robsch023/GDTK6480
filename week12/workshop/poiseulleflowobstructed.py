import numpy as np
import matplotlib.pyplot as plt
from time import time

DISCRETE_VELOCITIES = 9
DIRECTIONS = np.array([(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1),\
                       (1, 1), (-1, 1), (-1, -1), (1, -1)]).T
BOUNCED = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])

WEIGHTS = np.array([ 4/9, \
                     1/9,  1/9,  1/9,  1/9, \
                    1/36, 1/36, 1/36, 1/36], dtype=np.float64)
CS2 = 1/3  # lattice speed of sound squared


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

def f_equilibrium(rho, u_x, u_y):
    # f_i^eq calculation for the D2Q9 model
    f_eq = np.zeros((DISCRETE_VELOCITIES,) + rho.shape)
    for i in range(DISCRETE_VELOCITIES):
        cu = DIRECTIONS[0, i] * u_x + DIRECTIONS[1, i] * u_y
        f_eq[i] = rho * WEIGHTS[i] * (1 + cu/CS2 + cu**2/(2*CS2**2) - (u_x**2 + u_y**2)/(2*CS2))
    return f_eq


def moments(f, F_x=0, F_y=0):
    rho = np.sum(f, axis=0)
    u_x = (np.einsum("ijk, i -> jk", f, DIRECTIONS[0, :]) + 0.5*F_x) / rho
    u_y = (np.einsum("ijk, i -> jk", f, DIRECTIONS[1, :]) + 0.5*F_y) / rho
    return rho, u_x, u_y


def guo_source(u_x, u_y, F_x, F_y):
    # Guo et al. forcing term for the D2Q9 model
    S = np.zeros((DISCRETE_VELOCITIES,) + u_x.shape)
    for i in range(DISCRETE_VELOCITIES):
        c_x, c_y = DIRECTIONS[:, i]
        cu = c_x*u_x + c_y*u_y
        S[i] = WEIGHTS[i] * (((c_x - u_x)/CS2 + cu*c_x/CS2**2)* F_x
                        +    ((c_y - u_y)/CS2 + cu*c_y/CS2**2)* F_y)
    return S

"""
## VERIFICATION
rho = np.ones((1,1))
u_x = 0.05*np.ones((1,1))
u_y = 0.0*np.ones((1,1))
f_eq = f_equilibrium(rho, u_x, u_y)

print("f_eq at the node:", np.round(f_eq[:,0,0], 5))

rho_m, u_x_m, u_y_m = moments(f_eq)
print("rho from moments:", np.round(rho_m[0,0], 5))
print("u_x from moments:", np.round(u_x_m[0,0], 5))
print("u_y from moments:", u_y_m[0,0]) # little bit of rounding error
"""

# show that collision or relaxation process changes the distribution
# functions but preserves mass and momentum
NX = 300; NY = 80
N_ITER = 50_000

RE = 100
U0 = 0.05
cyl_cx = NX //4; cyl_cy = NY//2 + 2
cyl_r = NY//10

nu = U0*2*cyl_r/RE
tau = 3 * nu + 0.5
omega = 1/tau

# define geometry
xx, yy = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
solid = np.zeros((NX,NY), dtype=bool)
solid[:,0] = True
solid[:,-1] = True
solid[(xx - cyl_cx)**2 + (yy - cyl_cy)**2 < cyl_r**2] = True

## Initialisation
rho = np.ones((NX, NY))
u_x = U0 * np.zeros((NX, NY)); u_y = np.zeros((NX, NY))
u_x[solid] = 0.0
f_prec = f_equilibrium(rho, u_x, u_y)
mass0 = np.sum(f_prec)

fig, ax = plt.subplots(figsize=(6, 4))
# F_X = rho*G_X; F_Y = np.zeros((NX, NY))

for iter in range(N_ITER+1):
    # 1) MACROSCOPIC VARIABLES
    rho, u_x, u_y = moments(f_prec)

    # if iter % 20 == 0:
    #    print(f"iter {iter}: total mass change = {np.sum(rho) - mass0:.2e}")
    #    ax.plot(np.arange(NX) - NX//2, rho[:, NY//2] - 1, label = f"t  = {iter}")

    if iter % 200 == 0:
        plt.clf()
        plt.imshow(np.sqrt(u_x**2 + u_y**2).T, origin='lower', cmap='viridis')
        circle = plt.Circle((cyl_cx, cyl_cy), cyl_r)
        plt.gca().add_artist(circle)
        plt.pause(0.1)

    # 2) EQUILIBRIUM
    f_eq = f_equilibrium(rho, u_x, u_y)

    # 3) COLLISION
    f_solid = f_prec[:, solid]
    # f_post = f_prec - omega*(f_prec - f_eq) #+ (1 - 0.5*omega)*guo_source(u_x, u_y, G_X, 0)
    f_post = collisionMRT(f_prec, S2=1-omega)
    f_post[:, solid] = f_solid

    # 4) STREAMING
    for i in range(DISCRETE_VELOCITIES):
        f_prec[i, :, :] = np.roll(f_post[i, :, :], DIRECTIONS[:, i], axis=(0, 1))

    # 5) BOUNCE BACK
    f_tmp = f_prec.copy()
    for i in range(DISCRETE_VELOCITIES):
        f_prec[i, solid] = f_tmp[BOUNCED[i], solid]

    # 5b) wet node boundaries
    rho_w = 1/(1-U0) * (f_prec[0, 0, 1:-1] + f_prec[2, 0, 1:-1] + 2*f_prec[3, 0, 1:-1] + f_prec[4, 0, 1:-1] + 2*f_prec[6, 0, 1:-1] + 2*f_prec[7, 0, 1:-1])
    f_prec[1,0,1:-1] = f_prec[3,0,1:-1] + 0.666666666666667*rho_w*U0
    f_prec[5,0,1:-1] = -0.5*f_prec[2,0,1:-1] + 0.5*f_prec[4,0,1:-1] + f_prec[7,0,1:-1] + 0.166666666666667*rho_w*U0
    f_prec[8,0,1:-1] = 0.5*f_prec[2,0,1:-1] - 0.5*f_prec[4,0,1:-1] + f_prec[6,0,1:-1] + 0.166666666666667*rho_w*U0

    f_prec[3, -1, 1:-1] = f_prec[3, -2, 1:-1]
    f_prec[7, -1, 1:-1] = f_prec[7, -2, 1:-1]
    f_prec[6, 0, 1:-1] = f_prec[6, -2, 1:-1]

y = np.arange(1, NY-1)
H = NY - 2
y_tilde = y - 0.5
u_exact = G_X/(2*nu) * y_tilde * (H - y_tilde)
u_lbm = u_x[NX//2, 1:-1]
error = np.sqrt((np.sum((u_lbm - u_exact)**2) / np.sum(u_exact**2)))
print(f"tau = {TAU}, relative L2 error = {error:.2e}")

ax.plot(u_exact, y, 'k-', label='analytical')
ax.plot(u_lbm, y, 'ro', label=f'LBM, tau={TAU}')
plt.legend()
plt.show()

    