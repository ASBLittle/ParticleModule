""" Module defines some typical drag models for momentum exchange between the particle and fluid."""

import numpy

class model(object):
    def __init__(self, coefficient):
        self.coefficient = coefficient
    def __call__(self, fluid_velocity, particle_velocity, *args, **kwargs):
        return (self.coefficient(fluid_velocity, particle_velocity,
                                 *args, **kwargs)
                *(fluid_velocity-particle_velocity))

def stokes_drag_coefficient(fluid_velocity, particle_velocity,
                            diameter, rho, fluid_viscosity, **kwargs):
    del fluid_velocity
    del particle_velocity
    """ Return Stokes drag force for particle parameters"""
    del kwargs ## supress argument unused warning
    return rho*18./diameter**2*fluid_viscosity

def stokes_drag (fluid_velocity, particle_velocity,
                 diameter, rho, fluid_viscosity, **kwargs):
    """Calculate Stokes drag."""
    return stokes_drag_coefficient(fluid_velocity,
                                   particle_velocity,
                                   diameter, rho, fluid_viscosity,
                                   **kwargs
                                  )*(fluid_velocity-particle_velocity)

def turbulent_drag_coefficient(fluid_velocity, particle_velocity,
                               diameter, rho, rho_f, **kwargs):
    """ Return turbulent drag force for particle parameters"""
    del kwargs, rho ## supress argument unused warning
    fluid_velocity = numpy.array(fluid_velocity)
    particle_velocity = numpy.array(particle_velocity)
    delta = numpy.sqrt(sum((fluid_velocity-particle_velocity)**2))
    return rho_f*0.44*3.0/32.0/diameter*delta

def turbulent_drag(fluid_velocity, particle_velocity, diameter, rho, rho_f,
                   fluid_viscosity, **kwargs):
    """Calculate turbulent drag."""
    return turbulent_drag_coefficient(fluid_velocity,
                                      particle_velocity,
                                      diameter, rho, rho_f,
                                      **kwargs
                                     )*(fluid_velocity-particle_velocity)

def transitional_drag_coefficient(fluid_velocity, particle_velocity,
                                  diameter, rho_f, fluid_viscosity=1.0e-3,
                                  **kwargs):
    """ Return transitional drag force for particle parameters

    Transition to turbulence is assumed to occur when Re>1000 for a particle
    Reynolds number,

    Re=rho_f*du*diameter/ mu
    """
    del kwargs ## supress argument unused warning
    fluid_velocity = numpy.array(fluid_velocity)
    particle_velocity = numpy.array(particle_velocity)

    delta = numpy.sqrt(sum((fluid_velocity-particle_velocity)**2))

    reynolds_no = rho_f*delta*diameter/ fluid_viscosity

    if reynolds_no < 1.0e-8:
        return rho_f*delta*(fluid_velocity-particle_velocity)
    if reynolds_no < 1000.0:
        c_d = (24.0/reynolds_no)*(1.0+0.15*reynolds_no**0.687)
        return rho_f*c_d*3.0/32.0/diameter*delta
    #otherwise
    return rho_f*0.44*3.0/32.0/diameter*delta

def transitional_drag(fluid_velocity, particle_velocity, diameter, rho_f,
                   fluid_viscosity, **kwargs):
    return transitional_drag_coefficient(fluid_velocity,
                                         particle_velocity,
                                         diameter, rho_f,
                                         fluid_viscosity,
                                         **kwargs
                                        )*(fluid_velocity-particle_velocity)
