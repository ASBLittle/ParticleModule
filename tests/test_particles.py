import Particles

import vtk

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

