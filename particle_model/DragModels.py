""" Module defines some typical drag models for momentum exchange between the particle and fluid."""

import numpy

def stokes_drag(fluid_velocity, particle_velocity, diameter, fluid_viscosity, **kwargs):
    """ Return Stokes drag force for particle parameters"""
    del kwargs ## supress argument unused warning
    return 9./(2.*diameter)*fluid_viscosity*(fluid_velocity-particle_velocity)

def turbulent_drag(fluid_velocity, particle_velocity, diameter, **kwargs):
    """ Return turbulent drag force for particle parameters"""
    del kwargs ## supress argument unused warning
    fluid_velocity = numpy.array(fluid_velocity)
    particle_velocity = numpy.array(particle_velocity)
    delta = numpy.sqrt(sum((fluid_velocity-particle_velocity)**2))
    return 0.44*3.0/32.0/diameter*delta*(fluid_velocity-particle_velocity)

def transitional_drag(fluid_velocity, particle_velocity, diameter, rho_f=1.0e3,
                      fluid_viscosity=1.0e-3, **kwargs):
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
        return delta*(fluid_velocity-particle_velocity)
    elif reynolds_no < 1000.0:
        c_d = (24.0/reynolds_no)*(1.0+0.15*reynolds_no**0.687)
        return c_d*3.0/32.0/diameter*delta*(fluid_velocity-particle_velocity)
    else:
        return 0.44*3.0/32.0/diameter*delta*(fluid_velocity-particle_velocity)
