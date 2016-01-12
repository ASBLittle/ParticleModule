""" Base module containing classes used at multiple levels."""

import libspud
import numpy

from particle_model import Options
from particle_model import DragModels
from particle_model import Parallel

class ParticleBase(object):
    """ An easily picklable base class for checkpointing and parallel computation. """

    def __init__(self, pos, vel, time=0.0, delta_t=1.0, phash=None):

        self.pos = pos
        self.vel = vel
        self.time = time
        self.delta_t = delta_t
        self.id=Parallel.particle_id(phash)

    def update(self):
        """ Core method updating the particle."""
        raise NotImplementedError

    def exchange(self):
        """ A helper function for parallel coding"""
        raise NotImplementedError

    def __eq__(self, obj):
        return self.id == obj.id

    def __hash__(self):
        return self.id()

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
        self.drag = kwargs.get('drag', DragModels.transitional_drag)

    def __call__(self, key='diameter'):
        """Get attributes."""
        if key == 'diameter':
            return self.diameter
        elif key == 'rho':
            return self.rho
        else:
            return self.data_dict[key]

    def get_volume(self):
        """Return particle volume."""
        return 1.0/6.0*numpy.pi*self.diameter**3

    def get_mass(self):
        """Return particle mass."""
        return self.rho*self.volume

    def randomize(self):
        """Update particle parameters from the given distribution"""

        if self.distribution:
            new_particle = copy.deepcopy(self)
            new_particle.diameter = self.distribution(self.base_diameter)
        else:
            new_particle = self

        return new_particle


def get_parameters_from_options(options_file=None,**kwargs):

    if options_file:
        libspud.load_options(options_file)

    parameters = []

    options_base = '/embedded_models/particle_model/particle_classes'

    def get_option(pclass,key,default=None):

        result = default
        if libspud.have_option('/'.join((options_base,pclass,key))):
            result = libspud.get_option('/'.join((options_base,pclass,key)))
        return result

    for i in range(libspud.get_number_of_children(options_base)):
        key = libspud.get_child_name(options_base,i)
        name = get_option(key,'name')

        diameter = get_option(key,'diameter')
        distribution = get_option(key,'distribution')
        density = get_option(key,'density',default=2.5e3)
        
        parameters.append(PhysicalParticle(diameter=diameter,
                                           rho=density,
                                           distribution=distribution,
                                           material_name=name,
                                           **kwargs))

    return parameters


def get_parameters_from_reader(reader,**kwargs):

    options_base = '/embedded_models/particle_model/particle_classes'

    parameters = []

    def get_option(pclass,key,default=None):
        result = default
        if libspud.have_option('/'.join((options_base,pclass,key))):
            result = libspud.get_option('/'.join((options_base,pclass,key)))
        return result

    for i in range(libspud.get_number_of_children(options_base)):
        key = libspud.get_child_name(options_base,i)
        name = get_option(key,'name')

        diameter = get_option(key,'diameter')
        distribution = get_option(key,'distribution')
        density = get_option(key,'density',default=2.5e3)
        
        parameters.append(PhysicalParticle(diameter=diameter,
                                           rho=density,
                                           distribution=distribution,
                                           material_name=name,
                                           **kwargs))


    return parameters

                      
        
    
