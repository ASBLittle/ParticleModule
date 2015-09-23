### Module defines some typical drag models for momentum exchange between the particle and fluid.

def stokes_drag(u,v,d,mu,**kwargs):
    return 9./(2.*d)*mu*(u-v)

def test_stokes_drag():
    assert stokes_drag(u=1.0,v=0.0,d=1.0,mu=1.0)==9.0/2.

def turbulent_drag(u,v,d,**kwargs):
    return 0.44*3.0/32.0/d*numpy.sqrt(sum((u-v)**2))*(u-v)
