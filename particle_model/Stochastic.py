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

# ASBL added
def split_3d_gaussian_drift(kappa_xy, kappa_z, pos, vel, time, delta_t):
    """ A Langevin equation style random walk term, to apply as a callback."""

    # Choose direction (some double counting going on)
    theta = numpy.random.uniform(0, 2.0*numpy.pi)
    phi = numpy.random.uniform(0, 2.0*numpy.pi)
    # Calculate stochastic part of magnitude.
    r = numpy.random.standard_normal()
    # Calculate time varying part of drift term.
    # dr ~ N(0, sqrt(2*d*kappa*dt))
    # but we calculate dr/dt
    K_xy = numpy.sqrt(4.0*kappa_xy)/numpy.sqrt(delta_t)
    K_z = numpy.sqrt(2.0*kappa_z)/numpy.sqrt(delta_t)

    return numpy.array((K_xy * r * numpy.cos(theta) * numpy.cos(phi),
                        K_xy * r * numpy.sin(theta) * numpy.cos(phi),
                        K_z * r * numpy.sin(phi)) )


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
        """ Nonconstant valued drift term"""

        return generic_gaussian_drift[dim](kappa, pos, vel, time, delta_t)

    return func

# ASBL Added - for split xy and z kappa specification
def split_constant_gaussian_drift_term(kappa_xy, kappa_z, dim=3):
    """ Factory for a split constant valued drift term"""

    def func(pos, vel, time, delta_t):
        """ split constant valued drift term"""

        return split_3d_gaussian_drift(kappa_xy, kappa_z, pos, vel, time, delta_t)

    return func

def generic_1d_gaussian_impulse(kappa, pos, vel, time, delta_t):
    """ A Langevin equation style random walk term, to apply as a callback."""
    # Calculate stochastic part of magnitude.
    r = numpy.random.standard_normal()
    # Calculate time varying part of impulse term.
    # dr ~ N(0, sqrt(2*d*kappa*dt))
    # but we calculate dr/dt
    del_t = max(delta_t, 1.0e-8)
    K = numpy.sqrt(2.0*kappa(pos))/(del_t*numpy.sqrt(del_t))

    return numpy.array((r*K, 0.0, 0.0))

def generic_2d_gaussian_impulse(kappa, pos, vel, time, delta_t):
    """ A Langevin equation style random walk term, to apply as a callback."""

    # Choose direction (magnitue is plus or minus, so only need 0 to pi.
    theta = numpy.random.uniform(0, numpy.pi)
    # Calculate stochastic part of magnitude.
    r = numpy.random.standard_normal()
    # Calculate time varying part of impulse term.
    # dr ~ N(0, sqrt(2*d*kappa*dt))
    # but we calculate dr/dt
    del_t = max(delta_t, 1.0e-8)
    K = numpy.sqrt(4.0*kappa(pos))/(del_t*numpy.sqrt(del_t))

    return numpy.array((r*numpy.cos(theta), r*numpy.sin(theta), 0.0))*K

def generic_3d_gaussian_impulse(kappa, pos, vel, time, delta_t):
    """ A Langevin equation style random walk term, to apply as a callback."""

    # Choose direction (some double counting going on)
    theta = numpy.random.uniform(0, 2.0*numpy.pi)
    phi = numpy.random.uniform(0, 2.0*numpy.pi)
    # Calculate stochastic part of magnitude.
    r = numpy.random.standard_normal()
    # Calculate time varying part of impulse term.
    # dr ~ N(0, sqrt(2*d*kappa*dt))
    # but we calculate dr/dt
    del_t = max(delta_t, 1.0e-8)
    K = numpy.sqrt(6.0*kappa(pos))/(del_t*numpy.sqrt(del_t))

    return K*numpy.array((r*numpy.cos(theta)*numpy.cos(phi),
                          r*numpy.sin(theta)*numpy.cos(phi),
                          r*numpy.sin(phi)))

# ASBL Added - split kappas but fixed in space.
def split_3d_gaussian_impulse(kappa_xy, kappa_z, pos, vel, time, delta_t):
    """ A Langevin equation style random walk term, to apply as a callback."""

    # Choose direction (some double counting going on)
    theta = numpy.random.uniform(0, 2.0*numpy.pi)
    phi = numpy.random.uniform(0, 2.0*numpy.pi)
    # Calculate stochastic part of magnitude.
    r = numpy.random.standard_normal()
    # Calculate time varying part of impulse term.
    # dr ~ N(0, sqrt(2*d*kappa*dt))
    # but we calculate dr/dt
    del_t = max(delta_t, 1.0e-8)
    K_xy = numpy.sqrt(4.0*kappa_xy)/(del_t*numpy.sqrt(del_t))
    K_z = numpy.sqrt(2.0*kappa_z)/(del_t*numpy.sqrt(del_t))

    return numpy.array((K_xy * r * numpy.cos(theta) * numpy.cos(phi),
                        K_xy * r * numpy.sin(theta) * numpy.cos(phi),
                        K_z * r * numpy.sin(phi)) )

generic_gaussian_impulse = [None,
                            generic_1d_gaussian_impulse,
                            generic_2d_gaussian_impulse,
                            generic_3d_gaussian_impulse]


def constant_gaussian_impulse_term(kappa, dim=3):
    """ Factory for a constant valued impulse term"""

    def func(pos, vel, time, delta_t):
        """ Constant valued impulse term"""
        return generic_gaussian_impulse[dim](lambda pos: kappa, pos, vel, time, delta_t)

    return func

# ASBL Added
def split_constant_gaussian_impulse_term(kappa_xy, kappa_z, dim=3):
    """ Factory for a constant valued impulse term"""

    def func(pos, vel, time, delta_t):
        """ Constant valued impulse term"""
        return split_3d_gaussian_impulse(kappa_xy, kappa_z, pos, vel, time, delta_t)

    return func

def nonconstant_gaussian_impulse_term(kappa, dim=3):
    """ Factory for a nonconstant valued impulse term"""

    def func(pos, vel, time, delta_t):
        """ Nononstant valued impulse term"""

        return generic_gaussian_impulse[dim](kappa, pos, vel, time, delta_t)

    return func
