### Module defines some typical drag models for momentum exchange between the particle and fluid.

import numpy

def stokes_drag(u,v,d,mu,**kwargs):
    return 9./(2.*d)*mu*(u-v)

def turbulent_drag(u,v,d,**kwargs):
    U=numpy.array(u)
    V=numpy.array(v)
    return 0.44*3.0/32.0/d*numpy.sqrt(sum((U-V)**2))*(U-V)



def turbulent_drag(u,v,d,**kwargs):
    U=numpy.array(u)
    V=numpy.array(v)
    return 0.44*3.0/32.0/d*numpy.sqrt(sum((U-V)**2))*(U-V)


def transitional_drag(u,v,d,**kwargs):
    U=numpy.array(u)
    V=numpy.array(v)

    du=numpy.sqrt(sum((U-V)**2))

    Re=1.0e3*du*d/1.0e-3

    if Re<1000.0:
        return (24.0/Re)*(1.0+0.15*Re**0.687)*3.0/32.0/d*du*(U-V)
    else:
        return 0.44*3.0/32.0/d*du*(U-V)
