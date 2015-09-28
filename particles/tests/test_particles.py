import Particles
import IO
import Collision
import DragModels

import vtk
import numpy

bnd=IO.boundaryData('tests/data/rightward_boundary.vtu')

def tc(fname='rightward_0.vtu'):
    def fn(x):
        r1=vtk.vtkXMLUnstructuredGridReader()
        r1.SetFileName('tests/data/'+fname)
        r1.Update()

        l1=vtk.vtkCellLocator()
        l1.SetDataSet(r1.GetOutput())
        l1.BuildLocator()

        return [[None,None,r1.GetOutput(),l1],
                [None,None,r1.GetOutput(),l1]], 0.0

    return fn


def test_tests():
    assert 1


def test_basic_particle_initialization():
    from numpy import zeros

    p=zeros(3)
    v=zeros(3)

    pt=Particles.particle(p,v)

    assert all(pt.p==p) and all(pt.v==v)



def test_basic_particle_bucket_initialization(tmpdir):
    from numpy import zeros

    N=10

    p=zeros((N,3))
    v=zeros((N,3))

    p=Particles.particle_bucket(p,v)


def test_particle_bucket_step_do_nothing(tmpdir):
    from numpy import zeros

    bndc=IO.boundaryData('tests/data/boundary_circle.vtu')

    N=1

    p=zeros((N,3))
    v=zeros((N,3))
    u=zeros((N,3))
    gp=zeros((N,3))

    pb=Particles.particle_bucket(p,v,0.0,dt=0.5,U=u,GP=gp,
                                 base_name='tests/data/circle',
                                 filename=tmpdir.join('data.dat').strpath,
                                 boundary=bndc)

    pb.run(5.0)

    assert pb.t==5.0
    assert all(pb.particles[0].p==0.0)
    assert all(pb.particles[0].v==0.0)


def test_picker():

    P=Particles.particle(0,0,tc=tc())

    u,gp=P.picker((0.5,0.5,0.0),0.0)


    assert all(u==numpy.array((1.0,0,0)))



def test_step_constant_velocity():

    p0=numpy.array((0.5,0.5,0.0))

    v0=numpy.array((1.0,0.0,0.0))

    P=Particles.particle(p0,v0,dt=0.1,diameter=numpy.infty,tc=tc(),boundary=bnd)
    P.update()
    assert all(P.p==numpy.array((0.6,0.5,0)))
    assert P.t==0.1
    P.update()
    assert all(P.p==numpy.array((0.7,0.5,0)))


def test_step_spin_up_turbulent_drag():

    p0=numpy.array((0.1,0.5,0.0))
    v0=numpy.array((0.0,0.0,0.0))

    P=Particles.particle(p0,v0,dt=0.001,tc=tc(),boundary=bnd,drag=DragModels.turbulent_drag)
    P.update()
    assert all(abs(P.p-numpy.array((0.100345,0.5,0)))<1.e-8)
    assert P.t==0.001

def test_step_spin_up_transitional_drag():

    p0=numpy.array((0.1,0.5,0.0))
    v0=numpy.array((0.0,0.0,0.0))

    P=Particles.particle(p0,v0,dt=0.001,tc=tc(),boundary=bnd)
    P.update()
    assert all(abs(P.p-numpy.array((0.10373956,0.5,0)))<1.e-8)
    assert P.t==0.001


def test_step_head_on_collision():

    p0=numpy.array((0.9995,0.5,0.0))
    v0=numpy.array((1.0,0.0,0.0))

    P=Particles.particle(p0,v0,dt=0.001,diameter=numpy.infty,tc=tc(),boundary=bnd,e=1.0)
    P.update()
    assert all(abs(P.p-numpy.array((0.9995,0.5,0)))<1.0e-8)
    assert all(P.v==numpy.array((-1.,0,0)))
    assert P.t==0.001

    assert len(P.collisions)==1
    assert all(P.collisions[0].x==numpy.array((1.,0.5,0.)))
    assert P.collisions[0].time==0.0005
    assert all(P.collisions[0].v==numpy.array((1.,0.,0.)))
    assert P.collisions[0].angle==numpy.pi/2.0

def test_diagonal_collision():

    p0=numpy.array((0.9995,0.4995,0.0))
    v0=numpy.array((1.0,1.0,0.0))

    P=Particles.particle(p0,v0,dt=0.001,diameter=numpy.infty,tc=tc(),boundary=bnd,e=1.0)
    P.update()
    assert all(abs(P.p-numpy.array((0.9995,0.5005,0)))<1.0e-8)
    assert all(P.v==numpy.array((-1.,1.0,0)))
    assert P.t==0.001

    assert len(P.collisions)==1
    assert all(P.collisions[0].x-numpy.array((1.,0.5,0.))<1.0e-8)
    assert P.collisions[0].time-0.0005<1e-8
    assert all(P.collisions[0].v==numpy.array((1.,1.,0.)))
    assert P.collisions[0].angle-numpy.pi/4.0<1e-10


def test_gyre_collision():

    ### Regression test for Mclaury coefficient

    bndg=IO.boundaryData('tests/data/gyre_boundary.vtu')

    from math import pi

    p0=numpy.array((0.8,0.45,0.0))
    v0=numpy.array((2.0*pi,0.0,0.0))

    P=Particles.particle(p0,v0,dt=0.001,diameter=100.0e-4,tc=tc('gyre_0.vtu'),boundary=bndg,e=1.0)

    for i in range(100):
        P.update()

    assert P.p[0]<1.0
    assert P.p[1]<1.0
    assert P.p[0]>0.0
    assert P.p[1]>0.0
    

    assert len(P.collisions)==1
    assert P.collisions[0].x[0]==1.0
    assert abs(Collision.MclauryMassCoeff(P.collisions[0])-0.17749677523046933)<1.0e-8


def test_coefficient_of_restitution():

    p0=numpy.array((0.95,0.5,0.0))
    v0=numpy.array((1.0,0.0,0.0))

    print bnd.bnd

    P=Particles.particle(p0,v0,dt=0.1,diameter=numpy.infty,tc=tc(),boundary=bnd,e=0.5)
    P.update()
    assert all(abs(P.p-numpy.array((0.975,0.5,0)))<1.0e-8)
    assert all(P.v==numpy.array((-0.5,0,0)))
    assert P.t==0.1

    assert len(P.collisions)==1
    assert all(P.collisions[0].x==numpy.array((1.,0.5,0.)))
    assert P.collisions[0].time==0.05
    assert all(P.collisions[0].v==numpy.array((1.,0.,0.)))
    assert P.collisions[0].angle==numpy.pi/2.0
