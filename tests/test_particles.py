import Particles

import vtk
import numpy

rb=vtk.vtkXMLUnstructuredGridReader()
rb.SetFileName('tests/data/rightward_boundary.vtu')
rb.Update()
bndl=vtk.vtkCellLocator()
bndl.SetDataSet(rb.GetOutput())
bndl.BuildLocator()

def tc(x):
    r1=vtk.vtkXMLUnstructuredGridReader()
    r1.SetFileName('tests/data/rightward_0.vtu')
    r1.Update()

    l1=vtk.vtkCellLocator()
    l1.SetDataSet(r1.GetOutput())
    l1.BuildLocator()

    return [[None,None,r1.GetOutput(),l1],
            [None,None,r1.GetOutput(),l1]], 0.0


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

    reader=vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName('tests/data/boundary_circle.vtu')
    reader.Update()

    print reader.GetOutput()

    bfile=reader.GetOutput()

    bndl=vtk.vtkCellLocator()
    bndl.SetDataSet(bfile)
    bndl.BuildLocator()

    N=1

    p=zeros((N,3))
    v=zeros((N,3))
    u=zeros((N,3))
    gp=zeros((N,3))

    pb=Particles.particle_bucket(p,v,0.0,dt=0.5,U=u,GP=gp,
                                 base_name='tests/data/circle',
                                 filename=tmpdir.join('data.dat').strpath,
                                 bndl=bndl)

    pb.run(5.0)

    assert pb.t==5.0
    assert all(pb.particles[0].p==0.0)
    assert all(pb.particles[0].v==0.0)


def test_picker():
    
    def tc(x):

        r1=vtk.vtkXMLUnstructuredGridReader()
        r1.SetFileName('tests/data/rightward_1.vtu')
        r1.Update()

        l1=vtk.vtkCellLocator()
        l1.SetDataSet(r1.GetOutput())
        l1.BuildLocator()

        return [[None,None,r1.GetOutput(),l1],
                [None,None,r1.GetOutput(),l1]], 0.0

    P=Particles.particle(0,0,tc=tc)

    u,gp=P.picker((0.5,0.5,0.0),0.0)


    assert all(u==numpy.array((1.0,0,0)))



def test_step_constant_velocity():
    rb=vtk.vtkXMLUnstructuredGridReader()
    rb.SetFileName('tests/data/rightward_boundary.vtu')
    rb.Update()
    bndl=vtk.vtkCellLocator()
    bndl.SetDataSet(rb.GetOutput())
    bndl.BuildLocator()

    def tc(x):
        r1=vtk.vtkXMLUnstructuredGridReader()
        r1.SetFileName('tests/data/rightward_0.vtu')
        r1.Update()

        l1=vtk.vtkCellLocator()
        l1.SetDataSet(r1.GetOutput())
        l1.BuildLocator()

        return [[None,None,r1.GetOutput(),l1],
                [None,None,r1.GetOutput(),l1]], 0.0

    p0=numpy.array((0.5,0.5,0.0))

    v0=numpy.array((1.0,0.0,0.0))

    P=Particles.particle(p0,v0,dt=0.1,tc=tc,bndl=bndl,bnd=rb.GetOutput())
    P.update()
    assert all(P.p==numpy.array((0.6,0.5,0)))
    assert P.t==0.1
    P.update()
    assert all(P.p==numpy.array((0.7,0.5,0)))


def test_step_spin_up():

    p0=numpy.array((0.1,0.5,0.0))
    v0=numpy.array((0.0,0.0,0.0))

    P=Particles.particle(p0,v0,dt=0.001,tc=tc,bndl=bndl,bnd=rb.GetOutput())
    P.update()
    assert all(abs(P.p-numpy.array((0.100345,0.5,0)))<1.e-8)
    assert P.t==0.001


def test_step_head_on_collision():

    p0=numpy.array((0.9995,0.5,0.0))
    v0=numpy.array((1.0,0.0,0.0))

    P=Particles.particle(p0,v0,dt=0.001,d=numpy.infty,tc=tc,bndl=bndl,bnd=rb.GetOutput(),e=1.0)
    P.update()
    assert all(abs(P.p-numpy.array((0.9995,0.5,0)))<1.0e-8)
    assert all(P.v==numpy.array((-1.,0,0)))
    assert P.t==0.001

    assert len(P.collisions)==1
    assert all(P.collisions[0].x==numpy.array((1.,0.5,0.)))
    assert P.collisions[0].t==0.0005
    assert all(P.collisions[0].v==numpy.array((1.,0.,0.)))
    assert P.collisions[0].angle==numpy.pi/2.0

def test_diagonal_collision():

    p0=numpy.array((0.9995,0.4995,0.0))
    v0=numpy.array((1.0,1.0,0.0))

    P=Particles.particle(p0,v0,dt=0.001,d=numpy.infty,tc=tc,bndl=bndl,bnd=rb.GetOutput(),e=1.0)
    P.update()
    assert all(abs(P.p-numpy.array((0.9995,0.5005,0)))<1.0e-8)
    assert all(P.v==numpy.array((-1.,1.0,0)))
    assert P.t==0.001

    assert len(P.collisions)==1
    assert all(P.collisions[0].x-numpy.array((1.,0.5,0.))<1.0e-8)
    assert P.collisions[0].t-0.0005<1e-8
    assert all(P.collisions[0].v==numpy.array((1.,1.,0.)))
    assert P.collisions[0].angle-numpy.pi/4.0<1e-10
    


def test_coefficient_of_restitution():

    p0=numpy.array((0.95,0.5,0.0))
    v0=numpy.array((1.0,0.0,0.0))

    P=Particles.particle(p0,v0,dt=0.1,d=numpy.infty,tc=tc,bndl=bndl,bnd=rb.GetOutput(),e=0.5)
    P.update()
    assert all(abs(P.p-numpy.array((0.975,0.5,0)))<1.0e-8)
    assert all(P.v==numpy.array((-0.5,0,0)))
    assert P.t==0.1

    assert len(P.collisions)==1
    assert all(P.collisions[0].x==numpy.array((1.,0.5,0.)))
    assert P.collisions[0].t==0.05
    assert all(P.collisions[0].v==numpy.array((1.,0.,0.)))
    assert P.collisions[0].angle==numpy.pi/2.0
