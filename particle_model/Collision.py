""" Module contains routines to deal with calculating wear from
collision information."""

import copy
import numpy

from particle_model.Debug import logger
from particle_model.ParticleBase import ParticleBase
import vtk

class BadCollisionException(Exception):
    """ Exception to deal with bad collisions"""
    pass

class OutletException(Exception):
    """ Exception to deal with collisions"""
    def __init__(self, pos_i, vel_i, cell_index=None, angle=None,
                 normal=None, delta_t=None):
        particle = ParticleBase(pos_i, vel_i)
        Exception.__init__(self, particle, pos_i, cell_index, angle,
                           normal, delta_t)
#        self.args = CollisionInfo(particle, pos_i, cell_index, angle, normal)
        self.pos = pos_i
        self.vel = vel_i
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

    def get_wear(self, wm, ER):
        """ Calculate wear induced by this collision"""

        #print "Lily collision wm %s, ER %s, rho %s, volume %s" % (wm, ER, self.rho, self.volume)
        if wm == "Arabnejad":
            return arabnejad(self, ER)
        elif wm == "McLaury":
            return mclaury_mass_coeff(self, ER)
        else:
            print "Lily no wear option selected"   
            return arabnejad(self, ER)

    def write_data(self, poly_data):
        """ Write data to polyfile """
        del poly_data
        wear = self.get_wear(wm, ER)
        del wear
        raise NotImplementedError

STANDARD_MATERIAL = {'n': 2, 'k': 1., 'H':1., 'F_s': 1., 'F_B':1.}
CARBON_STEEL1018 = {'n': 2.41, 'C': 0.01, 'H':1.31e11, 'F_s': 0.5, 'K':0.2, 'den':7870., 'eps':4.2e11, 'U0':5., 'F_B':1}
STAINLESS_STEEL316 = {'n': 2.41, 'C': 0.0125, 'H':2.24e11, 'F_s': 0.5, 'K':0.2, 'den':8000., 'eps':1.4e11, 'U0': 5.7, 'F_B':1}
ALUMINIUM6061 = {'n': 2.41, 'C': 0.0045, 'H':0.31e11, 'F_s': 0.5, 'K':0.2, 'den':2700., 'eps':0.6e11, 'U0': 7.3, 'F_B':1}

def basic_mclaury_mass_coeff(collision, material=None):
    """ Wear rate coefficient of collision from Mclaury correlation"""
    material = material or STANDARD_MATERIAL

    n_exp = material['n']
    cutting_coeff = material['C']
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

    return cutting_coeff*hardness*sharpness_factor*penetration_factor*vel**n_exp*fun(collision.angle)


def mclaury_mass_coeff(collision, ER = None, material=None):
    """ Wear rate coefficient of collision from Mclaury correlation"""
    #material = material or STANDARD_MATERIAL
    #material = material or CARBON_STEEL1018
    material = material or STAINLESS_STEEL316
    #material = material or ALUMINIUM6061

    n_exp = material['n']
    cutting_coeff = material['C']
    hardness = material['H']
    sharpness_factor = material['F_s']
    penetration_factor = material['F_B']
    ratio_of_contact = material['K']
    threshold_vel = material['U0']
    deformation = material['eps']
    p_rho = collision.rho
    p_volume = collision.volume
 
    logger.info('collision angle: %s', collision.angle)
    logger.info('collision normal: %s', collision.normal)

    def fun(theta):
        """ Mclaury angle response function"""
        if numpy.tan(theta) > ratio_of_contact:
        #if numpy.tan(theta) > 1.0/3.0:
            return numpy.cos(theta)**2/3.0
        #otherwise
        return numpy.sin(2.0*theta)-3.0*numpy.sin(theta)**2

    vel = numpy.sqrt(numpy.dot(collision.vel, collision.vel))
    #vel0 = 0.1
    vel0 = threshold_vel
    #beta = 1.00
    beta = sharpness_factor/deformation

    ER_C = (p_rho*p_volume*cutting_coeff*sharpness_factor*penetration_factor*(vel**n_exp*fun(collision.angle)))/numpy.sqrt(hardness)
    #ER_D = max(0.0, beta*(vel*numpy.sin(collision.angle)-vel0)**2) always positive?
    ER_D = p_rho*p_volume*beta*max(0.0, (vel*numpy.sin(collision.angle)-vel0))**2

    if ER == "ER_C":
        #print "Lily in MERC: ", ER_C
        return ER_C
    elif ER == "ER_D":
        #print "Lily in MERD: ", ER_D
        return ER_D
    else:
        #print "Lily ER Mtotal: ", (ER_C+ER_D)
        return (ER_C + ER_D)

    """return (cutting_coeff*hardness*sharpness_factor*penetration_factor
            *(vel**n_exp*fun(collision.angle)
              +max(0.0, beta*(vel*numpy.sin(collision.angle)-vel0)**2)))"""

def arabnejad(collision, ER = None, material=None):
    """Wear rate coefficient of collision from Arabnejad literature"""
    
    #material = material or CARBON_STEEL1018
    material = material or STAINLESS_STEEL316
    #material = material or ALUMINIUM6061

    n_exp = material['n']
    cutting_coeff = material['C']
    hardness = material['H']
    sharpness_factor = material['F_s']
    ratio_of_contact = material['K']
    material_density = material['den']
    deformation = material['eps']
    threshold_vel = material['U0']
    p_rho = collision.rho
    p_volume = collision.volume
    
    logger.info('collision angle: %s', collision.angle)
    logger.info('collision normal: %s', collision.normal)

    def fun(theta):
        """Arabnejad angle response function"""
        if numpy.tan(theta) > ratio_of_contact:
            return numpy.cos(theta)**2
        #otherwise
        return (ratio_of_contact*numpy.sin(2.0*theta)-numpy.sin(theta)**2)/ratio_of_contact**2

    vel = numpy.sqrt(numpy.dot(collision.vel, collision.vel))
    vel0 = threshold_vel
    beta = (sharpness_factor)/(deformation)
    ER_C = (p_rho*p_volume*cutting_coeff*sharpness_factor*vel**n_exp*fun(collision.angle))/numpy.sqrt(hardness)
    ER_D = p_rho*p_volume*beta*max(0.0,(vel*numpy.sin(collision.angle)-vel0))**2

    if ER == "ER_C":
        #print "Lily in AERC: ", ER_C
        return ER_C
    elif ER == "ER_D":
        #print "Lily in AERD: ", ER_D
        return ER_D
    else:
        #print "Lily ER Atotal: ", (ER_C+ER_D)
        return (ER_C + ER_D)
        
    ##return (ER_C+ER_D)

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
