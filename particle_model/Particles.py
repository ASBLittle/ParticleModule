""" Baseline module for the package. Contains the main classes, particle and particle_bucket. """

from particle_model import IO
from particle_model import DragModels
from particle_model import Collision
from particle_model import System

import numpy
import vtk
import scipy.linalg as la
import itertools
import copy

ARGV = [0.0, 0.0, 0.0]
ARGI = vtk.mutable(0)
ARGR = vtk.mutable(0.0)

def invert(mat):
    """ Hard coded 2D matrix inverse."""
    if mat.shape == (2, 2):
        return (numpy.array(((mat[1, 1], -mat[0, 1]),
                             (-mat[1, 0], mat[0, 0])))
                /(mat[0, 0]*mat[1, 1]-mat[0, 1]*mat[1, 0]))
    else:
        return la.inv(mat)

class PhysicalParticle(object):
    """ Class describing the physical properties of a particle drawn from a known distribution."""

    def __init__(self, diameter=40.0e-6, rho=2.5e3,
                 distribution=None, material_name='Sand', **kwargs):
        """Initialize physical particle state."""

        self.diameter = diameter
        self.base_diameter = diameter
        self.rho = rho

        self.distribution = distribution
        self.material_name = material_name
        self.data_dict = kwargs
        self.drag = kwargs.get('drag', DragModels.transitional_drag)

    def __call__(self, key='diameter'):
        """Get attributes."""
        if key == 'diameter':
            return self.diameter
        elif key == 'rho':
            return self.rho
        else:
            return self.data_dict[key]

    def get_volume(self):
        """Return particle volume."""
        return 1.0/6.0*numpy.pi*self.diameter**3

    def get_mass(self):
        """Return particle mass."""
        return self.rho*self.get_volume()

    def randomize(self):
        """Update particle parameters from the given distribution"""

        if self.distribution:
            new_particle = copy.deepcopy(self)
            new_particle.diameter = self.distribution(self.base_diameter)
        else:
            new_particle = self

        return new_particle

class Particle(object):
    """Class representing a single Lagrangian particle with mass"""

    def __init__(self, pos, vel, time=0.0, delta_t=1.0,
                 parameters=PhysicalParticle(diameter=40e-6, rho=2.5e3),
                 system=System.System()):

        self.pos = pos
        self.vel = vel
        self.time = time
        self.delta_t = delta_t
        self.collisions = []
        self.parameters = parameters
        self.system = system

    def update(self):
        """Update the state of the particle to the next time level

        The method uses relatively simple RK4 time integration."""

        kap1 = (self.vel, self.force(self.pos,
                                     self.vel,
                                     self.time))

        step, col, vel = self.collide(kap1[0], 0.5 * self.delta_t,
                                      self.vel,
                                      force=kap1[1])

        kap2 = (self.vel + 0.5*self.delta_t * kap1[1],
                self.force(self.pos + step,
                           self.vel + vel,
                           self.time + 0.5 * self.delta_t))

        step, col, vel = self.collide(kap2[0], 0.5 * self.delta_t,
                                      vel=self.vel, force=kap2[1])
        kap3 = (self.vel+0.5 * self.delta_t * kap2[1],
                self.force(self.pos + step,
                           self.vel + vel,
                           self.time + 0.5*self.delta_t))

        step, col, vel = self.collide(kap3[0], self.delta_t,
                                      vel=self.vel,
                                      force=kap3[1])
        kap4 = (self.vel + self.delta_t * kap3[1],
                self.force(self.pos + step,
                           self.vel + vel,
                           self.time + self.delta_t))


        step, col, vel = self.collide((kap1[0] + 2.0 * (kap2[0] + kap3[0]) + kap4[0]) / 6.0,
                                      self.delta_t, vel=self.vel,
                                      force=(kap1[1] + 2.0 * (kap2[1] + kap3[1]) + kap4[1])/6.0)
        self.pos += step
        self.vel += vel
        if col:
            self.collisions += col


        self.time += self.delta_t


    def get_fluid_properties(self):
        """ Get the fluid velocity and pressure gradient at the particle
        location."""
        return self.picker(self.pos, self.time)

    def force(self, position, particle_velocity, time):
        """Calculate the sum of the forces on the particle.

        Args:
            p (float): Location at which forcing is evaluated.
            v (float): Particle velocity at which forcing is evaluated.
            t (float): Time at which particle is evaluated.
        """

        fluid_velocity, grad_p = self.picker(position, time)

#        if collision:
#            raise collisionException

        return (grad_p / self.parameters.rho
                + self.parameters.drag(fluid_velocity,
                                       particle_velocity,
                                       self.parameters.diameter,
                                       fluid_viscosity=self.system.viscosity)
                + self.coriolis_force(particle_velocity)
                + self.system.gravity
                + self.centrifugal_force(position))

    def coriolis_force(self, particle_velocity):
        """ Return Coriolis force on particle."""

        def cross(vec1, vec2):
            """Return cross product of 3-tuples x and y."""
            out = numpy.zeros(3)
            out[0] = vec1[1]*vec2[2]-vec1[2]*vec2[1]
            out[1] = vec1[2]*vec2[0]-vec1[0]*vec2[2]
            out[2] = vec1[0]*vec2[1]-vec1[1]*vec2[0]
            return out

        return -2.0 * cross(self.system.omega, particle_velocity)

    def centrifugal_force(self, position):
        """ Return centrifugal force on particle"""
        return - (numpy.dot(self.system.omega, position) * self.system.omega
                  - numpy.dot(self.system.omega, self.system.omega) * position)

    def find_cell(self, locator, point=None):
        """ Use vtk rountines to find cell/element containing the point."""
        if point is None:
            point = self.pos
        cell_index = vtk.mutable(0)


        locator.FindClosestPoint(point, ARGV, cell_index, ARGI, ARGR)

        return cell_index

    def picker(self, pos, time):
        """ Extract fluid velocity and pressure from .vtu files at correct time level"""

        def fpick(infile, locator):
            """ Extract fluid velocity and pressure from single .vtu file"""

            locator.BuildLocatorIfNeeded()

            cell_index = self.find_cell(locator, pos)
            cell = infile.GetCell(cell_index)
            linear_cell = IO.get_linear_cell(cell)
            pids = cell.GetPointIds()

            dim = linear_cell.GetNumberOfPoints()-1
            upos = numpy.zeros(linear_cell.GetNumberOfPoints())
            dummy_func = linear_cell.GetPoints().GetPoint
            args = [dummy_func(i+1)[:dim] for i in range(dim)]
            args.append(dummy_func(0)[:dim])
            args.append(upos)
            linear_cell.BarycentricCoords(pos[:dim], *args)

#           collision == Collision.testInCell(linear_cell, pos)

            data_u = infile.GetPointData().GetVectors('Velocity')
            data_p = infile.GetPointData().GetScalars('Pressure')


            shape_funs = numpy.zeros(cell.GetNumberOfPoints())
            deriv_funs = numpy.zeros(dim*cell.GetNumberOfPoints())
            cell.InterpolateFunctions(upos[:3], shape_funs)
            cell.InterpolateDerivs(upos[:3], deriv_funs)


            rhs = numpy.array([data_p.GetValue(cell.GetPointId(i+1))
                               -data_p.GetValue(cell.GetPointId(0))
                               for i in range(dim)])

            mat = numpy.zeros((dim, dim))

            pts = numpy.array([cell.GetPoints().GetPoint(i) for i in range(dim+1)])
            for i in range(dim):
                mat[i, :] = pts[i+1, :dim] - pts[0, :dim]

#            mat=la.inv(mat)
            mat = invert(mat)

            nvout = numpy.array([data_u.GetTuple(pids.GetId(i))
                                 for i in range(cell.GetNumberOfPoints())])
            out = numpy.dot(shape_funs, nvout)

            grad_p = numpy.zeros(3)
            grad_p[:dim] = numpy.dot(mat, rhs)

            return out, grad_p

        data, alpha = self.system.temporal_cache(time)

        vel0, grad_p0 = fpick(data[0][2], data[0][3])
        vel1, grad_p1 = fpick(data[1][2], data[1][3])

        return ((1.0-alpha) * vel0 + alpha * vel1,
                (1.0 - alpha) * grad_p0 + alpha * grad_p1)


    def collide(self, k, delta_t, vel=None, force=None, pa=None, level=0):
        """Collision detection routine.

        Args:
            k  (float): Displacement
            dt (float): Timestep
            v  (float, optional): velocity
            f  (float, optional): forcing
            pa (float, optional): starting position in subcycle
            level (int) count to control maximum depth
        """
        if pa  is None:
            pa = self.pos

        if level == 10:
            return k * delta_t, None, vel - self.vel

        pos = pa+delta_t*k

        s = vtk.mutable(-1.0)
        x = [0.0, 0.0, 0.0]
        cell_index = vtk.mutable(0)

        bndl = self.system.boundary.bndl

        intersect = bndl.IntersectWithLine(pa, pos,
                                           1.0e-6, s,
                                           x, ARGV, ARGI, cell_index)

        if intersect:
            data, _ = self.system.temporal_cache(self.time)
            assert IO.test_in_cell(data[0][2].GetCell(self.find_cell(data[0][3], pa)), pa)

            print 'collision', intersect, cell_index, s, x, pos, pa
            x = numpy.array(x)

            cell = self.system.boundary.bnd.GetCell(cell_index)

            normal = numpy.zeros(3)

            vec1 = (numpy.array(cell.GetPoints().GetPoint(1))
                    -numpy.array(cell.GetPoints().GetPoint(0)))

            if cell.GetCellType() == vtk.VTK_TRIANGLE:
                vec2 = (numpy.array(cell.GetPoints().GetPoint(2))
                        -numpy.array(cell.GetPoints().GetPoint(0)))
            else:
                vec2 = numpy.array(((pos-pa)[1]*vec1[2]-(pos-pa)[2]*vec1[1],
                                    (pos-pa)[2]*vec1[0]-(pos-pa)[0]*vec1[2],
                                    (pos-pa)[0]*vec1[1]-(pos-pa)[1]*vec1[0]))

            normal[0] = vec1[1]*vec2[2]-vec1[2]*vec2[1]
            normal[1] = vec1[2]*vec2[0]-vec1[0]*vec2[2]
            normal[2] = vec1[0]*vec2[1]-vec1[1]*vec2[0]


            if sum(normal**2) > 1.0e-32:
                normal = normal / numpy.sqrt(sum(normal**2))
            else:
                print normal
                print vec1, vec2
                raise Collision.CollisionException

            normal = normal * numpy.sign(numpy.dot(normal, (pos-pa)))

            coeff = self.system.coefficient_of_restitution(self, cell)

            pos = x + delta_t * (k - (1.0 + coeff) * normal * (numpy.dot(normal, k)))

            theta = abs(numpy.arcsin(numpy.dot(normal, (x-pa))
                                     / numpy.sqrt(numpy.dot(x - pa, x - pa))))

            coldat = []

            if any(vel):
                vels = vel + s * delta_t * force
            else:
                vels = self.vel + s * delta_t * force

            par_col = copy.copy(self)
            par_col.pos = x
            par_col.vel = vels
            par_col.time = self.time + s * delta_t

            coldat.append(Collision.CollisionInfo(par_col, cell_index,
                                                  theta, normal))
            vels += -(1.0 + coeff)* normal * numpy.dot(normal, vels)

            px, col, velo = self.collide(vels, (1 - s) * delta_t,
                                         vel=vels, force=force,
                                         pa=x + 1.0e-10 * vels,
                                         level=level + 1)
            pos = px + x + 1.0e-10 * vels

            if col:
                coldat += col

            return pos - pa, coldat, velo

        return pos - pa, None, vel + delta_t * force - self.vel

class ParticleBucket(object):
    """Class for a container for multiple Lagrangian particles."""

    def __init__(self, X, V, time=0, delta_t=1.0e-3, filename=None,
                 parameters=PhysicalParticle(),
                 system=System.System()):
        """Initialize the bucket

        Args:
            X (float): Initial particle positions.
            V (float): Initial velocities
        """

        self.system = system

        self.fluid_velocity = numpy.zeros(X.shape)
        self.grad_p = numpy.zeros(X.shape)

        self.particles = []
        for _, (dummy_pos, dummy_vel) in enumerate(zip(X, V)):
            self.particles.append(Particle(dummy_pos, dummy_vel, time, delta_t,
                                           system=self.system,
                                           parameters=parameters.randomize()))
            if self.system.temporal_cache:
                self.fluid_velocity[_, :], self.grad_p[_, :] = \
                    self.particles[-1].get_fluid_properties()
        self.time = time
        self.pos = X
        self.vel = V
        self.delta_t = delta_t
        if filename:
            self.outfile = open(filename, 'w')

    def update(self):
        """ Update all the particles in the bucket to the next time level."""
        self.system.temporal_cache.range(self.time, self.time + self.delta_t)
        for k, part in enumerate(self.particles):
            part.update()
            if self.system.temporal_cache:
                self.fluid_velocity[k, :], self.grad_p[k, :] = \
                    part.get_fluid_properties()
        self.time += self.delta_t

    def collisions(self):
        """Collect all collisions felt by particles in the bucket"""
        return [i for i in itertools.chain(*[p.collisions for p in self.particles])]


    def write(self):
        """Write timelevel data to file."""

        self.outfile.write('%f'%self.time)

        for pos in self.pos.ravel():
            self.outfile.write(' %f'%pos)

        for vel in self.vel.ravel():
            self.outfile.write(' %f'%vel)

        self.outfile.write('\n')

    def run(self, time, write=True):
        """Drive particles forward until a given time."""
        while self.time - time < -1.0e-15:
            self.update()
            if write:
                self.write()
        if write:
            self.outfile.flush()
