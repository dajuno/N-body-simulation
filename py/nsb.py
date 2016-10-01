import numpy as np


class Body:
    def __init__(self, pos, vel, mass, id=None):
        ''' Single body class.

        Args:
            pos     initial position [x, y(, z)]
            vel     initial velocity vector [u_x, u_y(, u_v)]
            mass    mass
        '''
        self.x = np.array(pos)
        self.u = np.array(pos)
        self.a = np.zeros(self.x.shape)

        self.id = id

        self.mass = mass
        return self


class NBS():
    def __init__(self, num_bodies, radius=1, ang_velocity=1, mass_mean=1,
                 dt=1.e-3, num_iter=100, soften=1.e-8, dim=2):
        ''' Main N-Body Simulation class.

        Args:
            num_bodies      number of bodies
            radius          max radius of bodies
            ang_velocity    angular velocity of the disk
            mass_mean       mean mass
            dt              time step
            num_iter        number of iterations
            soften          softening constant, regularization
            dim             space dimensions (2, 3)
        '''
        self.num_bodies = num_bodies
        self.dt = dt
        self.soften = soften
        self.num_iter = num_iter

        self.bodies = []
        self.build_bodies(num_bodies, radius, ang_velocity, mass_mean, dim)

        self.G = 1.  # gravity constant

        return self

    def build_bodies(self, num, rad_max, ang_velocity, mass_mean, dim):
        ''' Build initial distribution of bodies.

        Args:
            num             number of bodies
            rad_max         max radius
            ang_velocity    angular velocity
            dim             spatial dimension
        '''
        if dim not in (2, 3):
            raise Exception('Dimension must be 2 or 3!')
        assert dim == 2, '3D not yet implemented'

        angles = np.random.uniform(0, 2*np.pi, num)
        radii = np.random.uniform(0, rad_max, num)
        masses = np.random.rand(mass_mean, mass_mean*0.5, num)

        for ang, rad, mass in zip(angles, radii, masses):
            pos = [rad*np.sin(ang), rad*np.cos(ang)]
            vel = [rad*ang_velocity*np.cos(ang), -rad*ang_velocity*np.sin(ang)]
            self.bodies.append(Body(pos, vel, mass))

        return self

    def force(self, b1, b2):
        ''' Calculate the gravity force between two bodies.

        Args:
            b1      body 1, instance of Body class
            b2      body 2, instance of Body class

        Returns:
            force
        '''
        dist = b1.x - b2.x
        force = self.G*b1.m*b2.m*dist/(np.linalg.norm(dist)**3 + self.soften)
        return force





