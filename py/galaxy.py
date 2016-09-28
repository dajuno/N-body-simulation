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

x = np.random.rand(N, dim) - 0.5
m = mass_sun + np.random.rand(N, dim)
v = np.random.rand(N, dim) - 0.5

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

# x = np.array([[-0.2, 0.3], [0.3, -0.2]])
# v = np.array([[0.2, -0.1], [-0.3, 0.07]])
z = np.vstack((v, x))
dt = 1e-2

plt.figure()
plt.ion()
plt.hold(False)
for i in range(Niter):
    v = z[0:N, :]
    x = z[N:, :]
    F = forces(x, m, G)
    Fv = sumgravity(F)
    z = z + dt*np.vstack((Fv, v))

    plt.plot(x[:, 0], x[:, 1], 'o')
    plt.axis([-1, 1, -1, 1])
    plt.draw()



