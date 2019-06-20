"""__init__.py file for the particle_model package"""
import particle_model.Debug
import particle_model.Parallel
import particle_model.Collision
import particle_model.IO
import particle_model.Particles
import particle_model.Coupling
import particle_model.DragModels
import particle_model.Stochastic
from particle_model.vtk_extras import __version__

from numpy import zeros

def setup_from_fluidity(mblock, time, system_dt, delta_t=None,
                        positions=None, velocities=None):
    """Drive particles from inside Fluidity."""

    positions = positions or zeros((0, 3))
    velocities = velocities or zeros((len(positions), 3))

    system = particle_model.System.get_system_from_options(block=(mblock, time, system_dt))
    parameters = particle_model.ParticleBase.get_parameters_from_options()[0]

    bucket = particle_model.Particles.ParticleBucket(positions, velocities,
                                                     time, delta_t,
                                                     system=system,
                                                     parameters=parameters)

    return bucket


def DebugOn():
    """ Turn on debug level logging."""
    logging = particle_model.Debug.logging
    logger = particle_model.Debug.logger
    logger.setLevel(level=logging.DEBUG)

def InfoOn():
    """ Turn on info level logging."""
    logging = particle_model.Debug.logging
    logger = particle_model.Debug.logger
    logger.setLevel(level=logging.INFO)
