import numpy as np
import matplotlib.pyplot as plt


def forces(x, m, G):
    # NOTE: F is antisymmetric!
    N = x.shape[0]
    dim = x.shape[1]
    F = np.zeros((N, N, dim))
    for i in range(len(F)):
        for j in range(len(F)):
            r = x[j, :] - x[i, :]
            F[i, j, :] = m[i]*m[j]*r/np.linalg.norm(r)**3
    F *= G
    return F

def sumgravity(F):
    ''' force pair interaction matrix F is nan on diagonal -- sum all other
    elements in row to get corresponding resulting force
    '''
    # TODO: use antisymmetricity of F

    N = F.shape[0]
    dim = F.shape[2]
    Fv = np.zeros((N, dim))

    # replace nans by 0
    F[np.isnan(F)] = 0
    Fv = np.sum(F, 1)
    return Fv

G = 6.673e-3  # N m^2 kg^-2
mass_sun = 1   # 1988500e24
# F = G*m1*m2/r**2

dim = 2
Niter = 1000

# initialize stars
N = 5

# random
x = np.random.rand(N, dim) - 0.5
m = mass_sun + np.random.rand(N, dim)
v = np.random.rand(N, dim) - 0.5

# square/circle (N=4)
N = 4
x = np.array([[-0.5, 0],
              [0, 0.5],
              [0.5, 0],
              [0, -0.5]])

v = np.array([[0., 1],
              [1, 0],
              [0, -1],
              [-1, 0]])

m = 20*np.ones(N)

# disk, polar coordinates
N = 120
#  r = np.ones(N)
r = np.random.rand(N)*10
#  alpha = np.linspace(0, 2*np.pi, N)
alpha = np.random.rand(N)*2*np.pi
omega = 1.
x = 0.5*(r*np.array([np.cos(alpha), np.sin(alpha)])).T
v = (omega*r*np.array([-np.sin(alpha), np.cos(alpha)])).T
m = 5*np.ones(N)


x_ode = x.copy()
v_ode = v.copy()
m_ode = m.copy()


z_ode = np.vstack((v_ode, x_ode))
# x = np.array([[-0.2, 0.3], [0.3, -0.2]])
# v = np.array([[0.2, -0.1], [-0.3, 0.07]])
#  z = np.vstack((v, x))
dt = 1e-2

plt.figure()
plt.ion()
plt.hold(False)
t = dt
x_res = np.zeros((Niter, dim))
t_res = np.zeros(Niter)
for i in range(Niter):
    #  F = forces(x, m, G)
    #  Fv = sumgravity(F)
    #  a = Fv
    #  x += 0.5*dt**2 + v*dt
    #  v += a*dt

    Fv_ode = sumgravity(forces(x_ode, m_ode, G))
    v_ode = z_ode[0:N, :]
    x_ode = z_ode[N:, :]
    z_ode = z_ode + dt*np.vstack((Fv_ode, v_ode))

    t += dt
    t_res[i] = t

    #  plt.plot(x[:, 0], x[:, 1], 'o')
    plt.plot(x_ode[:, 0], x_ode[:, 1], 'ro')
    plt.axis([-10, 10, -10, 10])
    plt.draw()

    t += dt
