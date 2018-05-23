""" Module containing functions which implement stochastic behaviour, e.g. for dispersion."""

import numpy
import numpy.random

def generic_1d_gaussian_drift(kappa, pos, vel, time, delta_t):
    """ A Langevin equation style random walk term, to apply as a callback."""
    # Calculate stochastic part of magnitude.
    r = numpy.random.standard_normal()
    # Calculate time varying part of drift term.
    # dr ~ N(0, sqrt(2*d*kappa*dt))
    # but we calculate dr/dt
    K = numpy.sqrt(2.0*kappa(pos))/numpy.sqrt(delta_t)

    return numpy.array((r*K, 0.0, 0.0))

def generic_2d_gaussian_drift(kappa, pos, vel, time, delta_t):
    """ A Langevin equation style random walk term, to apply as a callback."""

    # Choose direction (magnitue is plus or minus, so only need 0 to pi.
    theta = numpy.random.uniform(0, numpy.pi)
    # Calculate stochastic part of magnitude.
    r = numpy.random.standard_normal()
    # Calculate time varying part of drift term.
    # dr ~ N(0, sqrt(2*d*kappa*dt))
    # but we calculate dr/dt
    K = numpy.sqrt(4.0*kappa(pos))/numpy.sqrt(delta_t)

    return numpy.array((r*numpy.cos(theta), r*numpy.sin(theta), 0.0))*K

def generic_3d_gaussian_drift(kappa, pos, vel, time, delta_t):
    """ A Langevin equation style random walk term, to apply as a callback."""

    # Choose direction (some double counting going on)
    theta = numpy.random.uniform(0, 2.0*numpy.pi)
    phi = numpy.random.uniform(0, 2.0*numpy.pi)
    # Calculate stochastic part of magnitude.
    r = numpy.random.standard_normal()
    # Calculate time varying part of drift term.
    # dr ~ N(0, sqrt(2*d*kappa*dt))
    # but we calculate dr/dt
    K = numpy.sqrt(6.0*kappa(pos))/numpy.sqrt(delta_t)

    return K*numpy.array((r*numpy.cos(theta)*numpy.cos(phi), 
                          r*numpy.sin(theta)*numpy.cos(phi),
                          r*numpy.sin(phi)))

generic_gaussian_drift = [None,
                          generic_1d_gaussian_drift,
                          generic_2d_gaussian_drift,
                          generic_3d_gaussian_drift]


def constant_gaussian_drift_term(kappa, dim=3):
    """ Factory for a constant valued drift term"""

    def func(pos, vel, time, delta_t):
        """ Constant valued drift term"""
        return generic_gaussian_drift[dim](lambda pos: kappa, pos, vel, time, delta_t)
    
    return func

def nonconstant_gaussian_drift_term(kappa, dim=3):
    """ Factory for a nonconstant valued drift term"""

    def func(pos, vel, time, delta_t):
        """ Nononstant valued drift term"""

        return generic_gaussian_drift[dim](kappa, pos, vel, time, delta_t)
    
    return func
