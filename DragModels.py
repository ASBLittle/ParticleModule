### Module defines some typical drag models for momentum exchange between the particle and fluid.

import numpy

def stokes_drag(u,v,d,mu,**kwargs):
    return 9./(2.*d)*mu*(u-v)

def turbulent_drag(u,v,d,**kwargs):
    U=numpy.array(u)
    V=numpy.array(v)
    return 0.44*3.0/32.0/d*numpy.sqrt(sum((U-V)**2))*(U-V)
