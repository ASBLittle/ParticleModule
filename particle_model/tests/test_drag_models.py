import particle_model.DragModels as DM
from numpy import array

def test_stokes_drag():
    assert DM.stokes_drag(u=1.0,v=0.0,d=1.0,mu=1.0)==9.0/2.


def test_turbulent_drag():
    u=array((1.,0.,0.))
    v=array((0.,0.,0.))
    e=array((1.0e-8,1.0e-8,1.0e-8))
    q=0.44*3.0/32.0*u
    assert all(DM.turbulent_drag(u=u,v=v,d=1.0,mu=1.0)-q<e)

    u=array((0.,5.,0.))
    v=array((0.,1.,0.))
    q=0.44*3.0/20.0*v

    assert all(DM.turbulent_drag(u=u,v=v,d=10.0,mu=1.0)-q<e)
