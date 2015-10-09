""" Module defines some typical drag models for momentum exchange between the particle and fluid."""

import numpy

def stokes_drag(u, v, d, mu, **kwargs):
    """ Return Stokes drag force for particle parameters"""
    del kwargs ## supress argument unused warning
    return 9./(2.*d)*mu*(u-v)

def turbulent_drag(u, v, d, **kwargs):
    """ Return turbulent drag force for particle parameters"""
    del kwargs ## supress argument unused warning
    U=numpy.array(u)
    V=numpy.array(v)
    return 0.44*3.0/32.0/d*numpy.sqrt(sum((U-V)**2))*(U-V)

def transitional_drag(fluid_velocity,particle_velocity, diameter, rho_f=1.0e3, mu=1.0e-3, **kwargs):
    """ Return transitional drag force for particle parameters

    Transition to turbulence is assumed to occur when Re>1000 for a particle Reynolds number, 

    Re=rho_f*du*diameter/ mu
    
    """
    del kwargs ## supress argument unused warning
    U=numpy.array(fluid_velocity)
    V=numpy.array(particle_velocity)

    du=numpy.sqrt(sum((U-V)**2))

    Re=rho_f*du*diameter/ mu

    if Re<1.0e-8:
        return du*(U-V)
    elif Re<1000.0:
        return (24.0/Re)*(1.0+0.15*Re**0.687)*3.0/32.0/diameter*du*(U-V)
    else:
        return 0.44*3.0/32.0/diameter*du*(U-V)
