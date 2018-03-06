""" Module contains routines to deal with calculating wear from
collision information."""

import copy
import numpy

from particle_model.Debug import logger
import vtk

class BadCollisionException(Exception):
    """ Exception to deal with bad collisions"""
    pass

class OutletException(Exception):
    """ Exception to deal with collisions"""
    def __init__(self, particle, pos_i, cell_index=None, angle=None,
                 normal=None, delta_t=None):
        Exception.__init__(self, particle, pos_i, cell_index, angle,
                           normal, delta_t)
        self.args = CollisionInfo(particle, pos_i, cell_index, angle, normal)
        self.delta_t = delta_t

class MappedBoundaryException(Exception):
    """ Exception to deal with collisions"""
    def __init__(self, *args):
        Exception.__init__(self, *args)

class CollisionException(Exception):
    """ Exception to deal with collisions"""
    def __init__(self, particle, pos_i, cell_index, delta_t):
        angle, normal = collision_angle(particle, particle.pos, pos_i,
                                        cell_index)
        Exception.__init__(self, particle, pos_i, cell_index, delta_t)
        self.args = particle, pos_i, cell_index, delta_t
        self.pos = pos_i
        self.info = CollisionInfo(particle, pos_i, cell_index, angle, normal)
        self.delta_t = delta_t
        self.vel = None

class CollisionInfo(object):
    """ Utility class for collision information """
    def __init__(self, particle, pos, cell, angle, normal):
        """ Initialise from particle collision information."""
        self.particle = copy.copy(particle)
        self.pos = copy.deepcopy(pos)
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
        del wear
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
            return numpy.cos(theta)**2/3.0
        #otherwise
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

    logger.info('collision angle: %s', collision.angle)
    logger.info('collision normal: %s', collision.normal)

    def fun(theta):
        """ Mclaury angle response function"""
        if numpy.tan(theta) > 1.0/3.0:
            return numpy.cos(theta)**2/3.0
        #otherwise
        return numpy.sin(2.0*theta)-3.0*numpy.sin(theta)**2

    vel = numpy.sqrt(numpy.dot(collision.vel, collision.vel))
    vel0 = 0.1
    beta = 1.0

    return (coeff*hardness*sharpness_factor*penetration_factor
            *(vel**n_exp*fun(collision.angle)
              +max(0.0, beta*(vel*numpy.sin(collision.angle)-vel0)**2)))


def collision_angle(particle, pos_0, pos_i, cell_index):
    """ Calulate the surface normal and angle of incidence of a collision event."""

    cell = particle.system.boundary.bnd.GetCell(cell_index)

    normal = numpy.zeros(3)

    veci = pos_i-pos_0

    vec1 = (numpy.array(cell.GetPoints().GetPoint(1))
            -numpy.array(cell.GetPoints().GetPoint(0)))

    if cell.GetCellType() == vtk.VTK_TRIANGLE:
        vec2 = (numpy.array(cell.GetPoints().GetPoint(2))
                -numpy.array(cell.GetPoints().GetPoint(0)))
    else:
        vec2 = numpy.array((veci[1]*vec1[2]-veci[2]*vec1[1],
                            veci[2]*vec1[0]-veci[0]*vec1[2],
                            veci[0]*vec1[1]-veci[1]*vec1[0]))

    normal[0] = vec1[1]*vec2[2]-vec1[2]*vec2[1]
    normal[1] = vec1[2]*vec2[0]-vec1[0]*vec2[2]
    normal[2] = vec1[0]*vec2[1]-vec1[1]*vec2[0]

    if sum(normal**2) > 1.0e-32:
        normal = normal / numpy.sqrt(sum(normal**2))
    else:
        logger.error(normal)
        logger.error("%s, %s", vec1, vec2)
        raise BadCollisionException

    normal = normal * numpy.sign(numpy.dot(normal, -veci))

    theta = abs(numpy.arcsin(numpy.dot(normal, veci)
                             / numpy.sqrt(numpy.dot(veci, veci))))

    return theta, normal

def rebound_velocity(particle, vel_i, grid_vel, normal, cell_index):
    """ Calculate velocity post collision."""

    coeff = particle.system.coefficient_of_restitution(particle, cell_index)
    vel_o = vel_i
    vel_o += -(1.0 + coeff)* normal * numpy.dot(normal, vel_i-grid_vel)

    return vel_o
