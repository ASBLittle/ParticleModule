"""Module describing the basic fluids system in which the particles are embedded."""

from numpy import zeros
from numpy.linalg import norm

class System(object):
    """ Class decribes the fixed properties of the underlying system and its
    fluid dynamical solution."""

    def __init__(self, boundary=None, **kwargs):
        """ Initialise the system class. """
        self.boundary = boundary

        self.gravity = kwargs.get('gravity', zeros(3))
        self.omega = kwargs.get('omega', zeros(3))
        self.rho = kwargs.get('rho', 1.0e3)
        self.viscosity = kwargs.get('viscosity', 1.0e-3)
        self.coeff = kwargs.get('coeff', 0.99)


    def get_rossby(self, velocity, length):
        """ Get the rossby number of the system for a given length scale
        and velocity."""
        return norm(velocity)/(length * 2.0 * norm(self.omega))


    def get_reynolds(self, velocity, length):
        """ Get the Reynolds number of the system for a given lenght scale
        and velocity."""
        return self.rho*norm(velocity)*length/ self.viscosity

    def coefficient_of_restitution(self, particle, cell=None):
        """ Get coefficent of restitution."""

        del particle
        del cell

        return self.coeff
