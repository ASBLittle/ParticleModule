""" Base module containing classes used at multiple levels."""

import copy
try:
    import libspud
except ImportError:
    pass
import numpy

from particle_model import DragModels
from particle_model import Parallel

class ParticleBase(object):
    """ An easily picklable base class for checkpointing and parallel computation. """

    def __init__(self, pos, vel, time=0.0, delta_t=1.0, phash=None):

        self.pos = pos
        self.vel = vel
        self.time = time
        self.delta_t = delta_t
        self._hash = Parallel.ParticleId(phash)
        self.fields = {}
        self._old = []

    def __hash__(self):
        return hash(self._hash)

    def set_hash(self, phash):
        """Update particle hash."""
        self._hash = phash

    def update(self, delta_t, method):
        """ Core method updating the particle."""
        raise NotImplementedError

    def __eq__(self, obj):
        return hash(self) == hash(obj)

    def set_old(self, old, num_time_levels=1):
        """Update old particle data."""
        if self._old:
            self._old = [None]+self._old[0:num_time_levels-1]
        else:
            self._old = [None]
        self._old[0] = copy.deepcopy(old)

    def get_old(self, time_level, key=None):
        """ Get old particle data"""
        if key is None:
            return self._old[time_level]
        #otherwise
        return self._old[time_level][key]

class PhysicalParticle(object):
    """ Class describing the physical properties of a particle drawn from a known distribution."""

    def __init__(self, diameter=40.0e-6, rho=2.5e3,
                 distribution=None, material_name='Sand', **kwargs):
        """Initialize physical particle state."""

        self.diameter = diameter
        self.base_diameter = diameter
        self.rho = rho

        self.distribution = distribution
        self.material_name = material_name
        self.data_dict = kwargs
        self.drag = kwargs.get('drag',
                               DragModels.Model(DragModels.transitional_drag_coefficient))
        self.drag_coefficient = kwargs.get('drag_coefficient',
                                           DragModels.transitional_drag_coefficient)

    def __call__(self, key='diameter'):
        """Get attributes."""
        if key == 'diameter':
            return self.diameter
        if key == 'rho':
            return self.rho
        #otherwise
        return self.data_dict[key]

    def pure_lagrangian(self):
        """Test if particle is purely Lagrangian (ie of zero particle diameter)"""
        return self.base_diameter == 0

    def get_area(self):
        """Return particle volume."""
        return 1.0/4.0*numpy.pi*self.diameter**2

    def get_volume(self):
        """Return particle volume."""
        return 1.0/6.0*numpy.pi*self.diameter**3

    def get_mass(self):
        """Return particle mass."""
        return self.rho*self.get_volume()

    def randomize(self):
        """Update particle parameters from the given distribution"""

        if self.distribution:
            new_particle = copy.deepcopy(self)
            new_particle.diameter = self.distribution(self.base_diameter)
        else:
            new_particle = self

        return new_particle


def get_parameters_from_options(options_file=None, **kwargs):
    """Read particle data from Fluidity options file."""

    if options_file:
        libspud.load_options(options_file)

    parameters = []

    options_base = '/embedded_models/particle_model/particle_classes'

    def get_option(pclass, key, default=None):
        """Get option from key."""

        result = default
        if libspud.have_option('/'.join((options_base, pclass, key))):
            result = libspud.get_option('/'.join((options_base, pclass, key)))
        return result

    for i in range(libspud.get_number_of_children(options_base)):
        key = libspud.get_child_name(options_base, i)
        name = get_option(key, 'name')

        diameter = get_option(key, 'diameter')
        distribution = get_option(key, 'distribution')
        density = get_option(key, 'density', default=2.5e3)

        parameters.append(PhysicalParticle(diameter=diameter,
                                           rho=density,
                                           distribution=distribution,
                                           material_name=name,
                                           **kwargs))

    return parameters


def get_parameters_from_reader(reader, **kwargs):
    """Get particle parameter data from an option reader object."""

    del reader

    options_base = '/embedded_models/particle_model/particle_classes'

    parameters = []

    def get_option(pclass, key, default=None):
        """libspud wrapper."""
        result = default
        if libspud.have_option('/'.join((options_base, pclass, key))):
            result = libspud.get_option('/'.join((options_base, pclass, key)))
        return result

    for i in range(libspud.get_number_of_children(options_base)):
        key = libspud.get_child_name(options_base, i)
        name = get_option(key, 'name')

        diameter = get_option(key, 'diameter')
        distribution = get_option(key, 'distribution')
        density = get_option(key, 'density', default=2.5e3)

        parameters.append(PhysicalParticle(diameter=diameter,
                                           rho=density,
                                           distribution=distribution,
                                           material_name=name,
                                           **kwargs))


    return parameters
