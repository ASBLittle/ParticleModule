### Module defines some typical drag models for momentum exchange between the particle and fluid.

def stokes_drag(u,v,d,**kwargs):
    return 0.44*3.0/32.0/d*numpy.sqrt(sum((u-v)**2))*(u-v)
