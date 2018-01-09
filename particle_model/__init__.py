"""__init__.py file for the particle_model package"""
import particle_model.Debug
import particle_model.Collision
import particle_model.IO
import particle_model.Particles
import particle_model.Coupling
import particle_model.DragModels

from numpy import zeros

def setup_from_fluidity(mb, time, dt, delta_t=None,
                        X=zeros((0, 3)),
                        V=zeros((0, 3))):
    """Drive particles from inside Fluidity."""

    SYSTEM = particle_model.System.get_system_from_options(block=(mb, time, dt))
    PAR = particle_model.ParticleBase.get_parameters_from_options()[0]

    bucket = particle_model.Particles.ParticleBucket(X, V, time, delta_t,
                                                     system=SYSTEM,
                                                     parameters=PAR)

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
