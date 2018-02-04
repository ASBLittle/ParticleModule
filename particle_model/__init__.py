"""__init__.py file for the particle_model package"""
import particle_model.Debug
import particle_model.Parallel
import particle_model.Collision
import particle_model.IO
import particle_model.Particles
import particle_model.Coupling
import particle_model.DragModels

from numpy import zeros

def setup_from_fluidity(mb, time, dt, delta_t=None,
                        positions=zeros((0, 3)),
                        velocities=zeros((0, 3))):
    """Drive particles from inside Fluidity."""

    system = particle_model.System.get_system_from_options(block=(mb, time, dt))
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
