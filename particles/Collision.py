""" Module contains routines to deal with calculating wear from
collision information."""

import numpy
import copy

class collisionException(Exception):
    """ Exception to deal with bad collisions"""
    pass

class collisionInfo(object):
    """ Utility class for collision information """
    def __init__(self, particle, cell, angle, time):
        """ Initialise from particle collision information."""
        self.particle = copy.copy(particle)
        self.x = copy.copy(particle.p)
        self.v = copy.copy(particle.v)
        self.cell = cell
        self.angle = angle
        self.time = copy.copy(time)

    def get_wear(self):
        """ Calculate wear induced by this collision"""
        return MclauryMassCoeff(self)

    def write_data(self, poly_data):
        """ Write data to polyfile """
        del poly_data
        wear = self.get_wear()
        print wear
        raise NotImplementedError

standard_material = {'n': 2, 'k': 1., 'H':1., 'F_s': 1., 'F_B':1.}

def testInCell(cell, position):
    """ Check if point is in vtk cell"""
    return cell.GetParametricDistance(position) == 0

def MclauryMassCoeff(collision, material=None):
    """ Wear rate coefficient of collision from Mclaury correlation"""
    material = material or standard_material

    n_exp = material['n']
    k = material['k']
    H = material['H']
    F_s = material['F_s']
    F_B = material['F_s']

    def fun(theta):
        """ Mclaury angle response function"""
        if numpy.tan(theta) > 1.0/3.0:
            return numpy.cos(theta)**2
        else:
            return numpy.sin(2.0*theta)-3.0*numpy.sin(theta)**2

    vel = numpy.sqrt(numpy.sum(collision.v**2))

    return k*H*F_s*F_B*vel**n_exp*fun(collision.angle)
