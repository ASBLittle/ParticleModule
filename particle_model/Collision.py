""" Module contains routines to deal with calculating wear from
collision information."""

import numpy
import copy

class CollisionException(Exception):
    """ Exception to deal with bad collisions"""
    pass

class CollisionInfo(object):
    """ Utility class for collision information """
    def __init__(self, particle, cell, angle, normal):
        """ Initialise from particle collision information."""
        self.particle = copy.copy(particle)
        self.pos = copy.deepcopy(particle.pos)
        self.vel = copy.deepcopy(particle.vel)
        self.time = copy.deepcopy(particle.time)
        self.cell = cell
        self.angle = copy.deepcopy(angle)
        self.normal = copy.deepcopy(normal)

    def get_wear(self):
        """ Calculate wear induced by this collision"""
        return mclaury_mass_coeff(self)

    def write_data(self, poly_data):
        """ Write data to polyfile """
        del poly_data
        wear = self.get_wear()
        print wear
        raise NotImplementedError

STANDARD_MATERIAL = {'n': 2, 'k': 1., 'H':1., 'F_s': 1., 'F_B':1.}

def basic_mclaury_mass_coeff(collision, material=None):
    """ Wear rate coefficient of collision from Mclaury correlation"""
    material = material or STANDARD_MATERIAL

    n_exp = material['n']
    coeff = material['k']
    hardness = material['H']
    sharpness_factor = material['F_s']
    penetration_factor = material['F_B']

    def fun(theta):
        """ Mclaury angle response function"""
        if numpy.tan(theta) > 1.0/3.0:
            return numpy.cos(theta)**2
        else:
            return numpy.sin(2.0*theta)-3.0*numpy.sin(theta)**2

    vel = numpy.sqrt(numpy.sum(collision.vel**2))

    return coeff*hardness*sharpness_factor*penetration_factor*vel**n_exp*fun(collision.angle)


def mclaury_mass_coeff(collision, material=None):
    """ Wear rate coefficient of collision from Mclaury correlation"""
    material = material or STANDARD_MATERIAL

    n_exp = material['n']
    coeff = material['k']
    hardness = material['H']
    sharpness_factor = material['F_s']
    penetration_factor = material['F_B']

    def fun(theta):
        """ Mclaury angle response function"""
        if numpy.tan(theta) > 1.0/3.0:
            return numpy.cos(theta)**2
        else:
            return numpy.sin(2.0*theta)-3.0*numpy.sin(theta)**2


    vel = numpy.sqrt(numpy.sum(collision.vel**2))

    vel0 = 0.0

    return coeff*hardness*sharpness_factor*penetration_factor*(vel**n_exp*fun(collision.angle)+max(0.0,vel*numpy.sin(collision.angle)-vel0))
