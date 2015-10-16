""" Test the main particle routines."""
from particle_model import Particles
from particle_model import IO
from particle_model import Collision
from particle_model import DragModels

import vtk
import numpy
from numpy import pi, sin, cos

BOUNDARY = IO.BoundaryData('particle_model/tests/data/rightward_boundary.vtu')

MESH = IO.GmshMesh()
MESH.read('particle_model/tests/data/Structured.msh')

def temp_cache(fname='rightward_0.vtu',ldir='particle_model/tests/data'):
    """Mock temporal cache."""
    def fun(time):
        """Factory function to mock a temporal cache."""
        del time
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(ldir+'/'+fname)
        reader.Update()

        locator = vtk.vtkCellLocator()
        locator.SetDataSet(reader.GetOutput())
        locator.BuildLocator()

        return [[0.0, fname, reader.GetOutput(), locator],
                [1.0, fname, reader.GetOutput(), locator]], 0.0

    return fun

def test_tests():
    """ Test test structure with a minimal test."""
    assert 1

def test_basic_particle_initialization():
    """Test basic particle initialization"""
    from numpy import zeros

    pres = zeros(3)
    vel = zeros(3)

    part = Particles.Particle(pres, vel)

    assert all(part.p == pres) and all(part.v == vel)



def test_basic_particle_bucket_initialization():
    """ Test initializing a particle bucket."""
    from numpy import zeros

    num = 10

    pres = zeros((num, 3))
    vel = zeros((num, 3))

    part = Particles.ParticleBucket(pres, vel)

    assert part


def test_particle_bucket_step_do_nothing(tmpdir):
    """ Test initializing a full particle bucket."""
    from numpy import zeros

    bndc = IO.BoundaryData('particle_model/tests/data/boundary_circle.vtu')

    num = 1

    pres = zeros((num, 3))
    vel = zeros((num, 3))
    fluid_vel = zeros((num, 3))
    grad_p = zeros((num, 3))

    bucket = Particles.ParticleBucket(pres, vel, 0.0, dt=0.5, U=fluid_vel, GP=grad_p,
                                      base_name='particle_model/tests/data/circle',
                                      filename=tmpdir.join('data.dat').strpath,
                                      boundary=bndc)

    bucket.run(5.0)

    assert bucket.time == 5.0
    assert all(bucket.particles[0].p == 0.0)
    assert all(bucket.particles[0].v == 0.0)


def test_picker_constant():
    """Test vtk picker."""

    part = Particles.Particle(0, 0, tc=temp_cache())
    fluid_velocity, grad_p = part.picker((0.5, 0.5, 0.0), 0.0)

    assert all(fluid_velocity == numpy.array((1.0, 0.0, 0.0)))
    assert all(grad_p == numpy.array((0.0, 0.0, 0.0)))


def test_picker_linear(tmpdir):
    """Test vtk picker."""

    X = ((0.5, 0.5, 0.0),
         (0.25,0.75,0.0))

    ERR = numpy.array((1.0e-8, 1.0e-8, 1.0e-8))
    FNAME = tmpdir.join('linear.vtu').strpath

    print FNAME

    def vel(x):
        return numpy.array(( x[0],
                              x[1], 0))

    def pres(x):
        return x[0]

    IO.make_unstructured_grid(MESH,vel,pres,0.0,FNAME)

    part = Particles.Particle(0, 0, tc=temp_cache('linear.vtu',
                                                  tmpdir.strpath))
    for point in X:

        fluid_velocity, grad_p = part.picker(point, 0.0)
        
        assert all(abs(fluid_velocity - vel(point)) < ERR)
        assert all(grad_p == numpy.array((1.0, 0.0, 0.0)))



def test_step_constant_velocity():
    """Test single step at constant velocity."""

    pos = numpy.array((0.5, 0.5, 0.0))
    vel = numpy.array((1.0, 0.0, 0.0))

    part = Particles.Particle(pos, vel, dt=0.1, diameter=numpy.infty,
                              tc=temp_cache(), boundary=BOUNDARY)
    part.update()
    assert all(part.p == numpy.array((0.6, 0.5, 0.0)))
    assert part.t == 0.1
    part.update()
    assert all(part.p == numpy.array((0.7, 0.5, 0.0)))


def test_step_spin_up_turbulent_drag():
    """Test turbulent drag function"""

    pos = numpy.array((0.1, 0.5, 0.0))
    vel = numpy.array((0.0, 0.0, 0.0))

    part = Particles.Particle(pos, vel, dt=0.001, tc=temp_cache(), boundary=BOUNDARY,
                              drag=DragModels.turbulent_drag)
    part.update()
    assert all(abs(part.p - numpy.array((0.100345, 0.5, 0))) < 1.e-8)
    assert part.t == 0.001

def test_step_spin_up_transitional_drag():
    """ Test transitional drag function."""

    pos = numpy.array((0.1, 0.5, 0.0))
    vel = numpy.array((0.0, 0.0, 0.0))

    part = Particles.Particle(pos, vel, dt=0.001, tc=temp_cache(), boundary=BOUNDARY)
    part.update()
    assert all(abs(part.p - numpy.array((0.10373956, 0.5, 0))) < 1.e-8)
    assert part.t == 0.001


def test_step_head_on_collision():
    """ Test a head-on collision."""

    pos = numpy.array((0.9995, 0.5, 0.0))
    vel = numpy.array((1.0, 0.0, 0.0))

    part = Particles.Particle(pos, vel, dt=0.001, diameter=numpy.infty,
                              tc=temp_cache(), boundary=BOUNDARY, e=1.0)
    part.update()
    assert all(abs(part.p - numpy.array((0.9995, 0.5, 0.0))) < 1.0e-8)
    assert all(part.v == numpy.array((-1., 0., 0.)))
    assert part.t == 0.001

    assert len(part.collisions) == 1
    assert all(part.collisions[0].pos == numpy.array((1., 0.5, 0.)))
    assert part.collisions[0].time == 0.0005
    assert all(part.collisions[0].vel == numpy.array((1., 0., 0.)))
    assert part.collisions[0].angle == numpy.pi/2.0

def test_diagonal_collision():
    """Test a collision at an angle"""

    pos = numpy.array((0.9995, 0.4995, 0.0))
    vel = numpy.array((1.0, 1.0, 0.0))

    part = Particles.Particle(pos, vel, dt=0.001, diameter=numpy.infty,
                              tc=temp_cache(), boundary=BOUNDARY, e=1.0)
    part.update()
    assert all(abs(part.p - numpy.array((0.9995, 0.5005, 0))) < 1.0e-8)
    assert all(part.v == numpy.array((-1., 1.0, 0.0)))
    assert part.t == 0.001

    assert len(part.collisions) == 1
    assert all(part.collisions[0].pos - numpy.array((1., 0.5, 0.)) < 1.0e-8)
    assert part.collisions[0].time - 0.0005 < 1e-8
    assert all(part.collisions[0].vel == numpy.array((1., 1., 0.)))
    assert part.collisions[0].angle - numpy.pi / 4.0 < 1e-10


def test_gyre_collision():
    """Regression test for Mclaury coefficient"""

    bndg = IO.BoundaryData('particle_model/tests/data/gyre_boundary.vtu')

    from math import pi

    pos = numpy.array((0.8, 0.45, 0.0))
    vel = numpy.array((2.0 * pi, 0.0, 0.0))

    part = Particles.Particle(pos, vel, dt=0.001, diameter=100.0e-4,
                              tc=temp_cache('gyre_0.vtu'), boundary=bndg, e=1.0)

    for i in range(100):
        del i
        part.update()

    assert part.p[0] < 1.0
    assert part.p[1] < 1.0
    assert part.p[0] > 0.0
    assert part.p[1] > 0.0

    assert len(part.collisions) == 1
    assert part.collisions[0].pos[0] == 1.0
    assert abs(Collision.mclaury_mass_coeff(part.collisions[0]) - 0.18205645627433897 ) < 1.0e-8


def test_coefficient_of_restitution():
    """Test of coefficient of restitution parameter."""

    pos = numpy.array((0.95, 0.5, 0.0))
    vel = numpy.array((1.0, 0.0, 0.0))

    part = Particles.Particle(pos, vel, dt=0.1, diameter=numpy.infty, tc=temp_cache(),
                              boundary=BOUNDARY, e=0.5)
    part.update()
    assert all(abs(part.p-numpy.array((0.975, 0.5, 0))) < 1.0e-8)
    assert all(part.v == numpy.array((-0.5, 0, 0)))
    assert part.t == 0.1

    assert len(part.collisions) == 1
    assert all(part.collisions[0].pos == numpy.array((1., 0.5, 0.)))
    assert part.collisions[0].time == 0.05
    assert all(part.collisions[0].vel == numpy.array((1., 0., 0.)))
    assert part.collisions[0].angle == numpy.pi/2.0
