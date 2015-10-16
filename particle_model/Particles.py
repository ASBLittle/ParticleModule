""" Baseline module for the package. Contains the main classes, particle and particle_bucket. """

from particle_model import TemporalCache
from particle_model import IO
from particle_model import DragModels
from particle_model import Collision

import numpy
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import vtk
import scipy.linalg as la
import itertools
import copy

ARGV = [0.0, 0.0, 0.0]
ARGI = vtk.mutable(0)
ARGR = vtk.mutable(0.0)

def invert(mat):
    if mat.shape == (2, 2):
        return numpy.array(((mat[1,1], -mat[0,1]), 
                            (-mat[1,0], mat[0,0])))/(mat[0,0]*mat[1,1]-mat[0,1]*mat[1,0])
    else:
        return la.inv(mat)

class Particle(object):
    """Class representing a single Lagrangian particle with mass"""

    def __init__(self, p, v, t=0.0, dt=1.0, tc=None, u=numpy.zeros(3),
                 gp=numpy.zeros(3), rho=2.5e3, g=numpy.zeros(3),
                 omega=numpy.zeros(3), diameter=40e-6, boundary=None,
                 e=0.99, drag=DragModels.transitional_drag):

        self.p = p
        self.v = v
        self.t = t
        self.dt = dt
        self.collisions = []
        self.tc = tc
        self.rho = 2.0e3
        self.g = g
        self.omega = omega
        self.diameter = diameter
        self.u = u
        self.gp = gp
        self.boundary = boundary
        self.e = e
        self.rho = rho
        self.drag = drag

    def update(self):
        """Update the state of the particle to the next time level

        The method uses relatively simple RK4 time integration."""

        kap1 = (self.v, self.force(self.p,
                                   self.v,
                                   self.t))

        step, col, vel = self.collide(kap1[0], 0.5 * self.dt, self.v, f=kap1[1])
        kap2 = (self.v + 0.5*self.dt * kap1[1], self.force(self.p + step,
                                                           self.v + vel,
                                                           self.t + 0.5 * self.dt))

        step, col, vel = self.collide(kap2[0], 0.5 * self.dt, v=self.v, f=kap2[1])
        kap3 = (self.v+0.5 * self.dt * kap2[1], self.force(self.p + step,
                                                           self.v + vel,
                                                           self.t + 0.5*self.dt))   

        step, col, vel = self.collide(kap3[0], self.dt, v=self.v, f=kap3[1])
        kap4 = (self.v + self.dt * kap3[1], self.force(self.p + step,
                                                       self.v + vel,
                                                       self.t + self.dt))


        step, col, vel = self.collide((kap1[0] + 2.0 * (kap2[0] + kap3[0]) + kap4[0]) / 6.0,
                                      self.dt, v=self.v,
                                      f=(kap1[1] + 2.0 * (kap2[1] + kap3[1]) + kap4[1])/6.0)
        self.p += step
        self.v += vel
        if col:
            self.collisions += col


        self.t += self.dt
        self.u[:], self.gp[:] = self.picker(self.p, self.t)

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

        return (grad_p / self.rho
                + self.drag(fluid_velocity, particle_velocity, self.diameter)
                + self.coriolis_force(particle_velocity)
                + self.g
                - self.centrifugal_force(position))

    def coriolis_force(self, particle_velocity):
        """ Return Coriolis force on particle."""

        def cross(x,y):
            """Return cross product of 3-tuples x and y."""
            z = numpy.zeros(3)
            z[0] = x[1]*y[2]-x[2]*y[1]
            z[1] = x[2]*y[0]-x[0]*y[2]
            z[2] = x[0]*y[1]-x[1]*y[0]
            return z

        return -2.0 * cross(self.omega, particle_velocity)

    def centrifugal_force(self, position):
        """ Return centrifugal force on particle"""
        return - ( numpy.dot(self.omega, position) * self.omega 
                   - numpy.dot(self.omega, self.omega) * position)

    def find_cell(self, locator, point=None):
        """ Use vtk rountines to find cell/element containing the point."""
        if point is None:
            point = self.p
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

            N = linear_cell.GetNumberOfPoints()-1
            x = numpy.zeros(linear_cell.GetNumberOfPoints())
            dummy_func=linear_cell.GetPoints().GetPoint
            args = [dummy_func(i)[:N] for i in range(1,N+1)]
            args.append(dummy_func(0)[:N])
            args.append(x)
            linear_cell.BarycentricCoords(pos[:N], *args)

#           collision == Collision.testInCell(linear_cell, pos)

            data_u = infile.GetPointData().GetVectors('Velocity')
            data_p = infile.GetPointData().GetScalars('Pressure')


            sf = numpy.zeros(cell.GetNumberOfPoints())
            df = numpy.zeros(N*cell.GetNumberOfPoints())
            cell.InterpolateFunctions(x[:3], sf)
            cell.InterpolateDerivs(x[:3], df)


            rhs = numpy.array([data_p.GetValue(cell.GetPointId(i+1)) - data_p.GetValue(cell.GetPointId(0)) for i in range(N)])

            mat = numpy.zeros((N, N))

            pts = numpy.array([cell.GetPoints().GetPoint(i) for i in range(N+1)])
            for i in range(N):
                mat[i, :] = pts[i+1,:N] - pts[0,:N]

#            mat=la.inv(mat)
            mat = invert(mat)

            nvout=numpy.array([data_u.GetTuple(pids.GetId(i)) 
                               for i in range(cell.GetNumberOfPoints())])
            out = numpy.dot(sf, nvout)

            grad_p = numpy.zeros(3)
            grad_p[:N] = numpy.dot(mat, rhs)

            return out, grad_p

        data, alpha = self.tc(time)

        vel0, grad_p0 = fpick(data[0][2], data[0][3])
        vel1, grad_p1 = fpick(data[1][2], data[1][3])

        return ((1.0-alpha) * vel0 + alpha * vel1,
                (1.0 - alpha) * grad_p0 + alpha * grad_p1)


    def collide(self, k, dt, v=None, f=None, pa=None, level=0):
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
            pa = self.p

        if level == 10:
            return k * dt, None, v - self.v

        p = pa+dt*k

        s = vtk.mutable(-1.0)
        x = [0.0, 0.0, 0.0]
        arg6 = [0.0, 0.0, 0.0]
        arg7 = vtk.mutable(0)
        cell_index = vtk.mutable(0)

        intersect = self.boundary.bndl.IntersectWithLine(pa, p,
                                                         1.0e-6, s,
                                                         x, arg6, arg7, cell_index)

        if s != -1.0:
            print 'collision', intersect,cell_index, s, x
            x = numpy.array(x)

            cell = self.boundary.bnd.GetCell(cell_index)

            normal = numpy.zeros(3)

            vec1 = (numpy.array(cell.GetPoints().GetPoint(1))
                    -numpy.array(cell.GetPoints().GetPoint(0)))

            if cell.GetCellType() == vtk.VTK_TRIANGLE:
                vec2 = (numpy.array(cell.GetPoints().GetPoint(2))
                        -numpy.array(cell.GetPoints().GetPoint(0)))
            else:
                
                vec2  = numpy.array((( x - pa )[1]*vec1[2]-( x - pa )[2]*vec1[1],
                                     ( x - pa )[2]*vec1[0]-( x - pa )[0]*vec1[2],
                                     ( x - pa )[0]*vec1[1]-( x - pa )[1]*vec1[0]))

            normal[0] = vec1[1]*vec2[2]-vec1[2]*vec2[1]
            normal[1] = vec1[2]*vec2[0]-vec1[0]*vec2[2]
            normal[2] = vec1[0]*vec2[1]-vec1[1]*vec2[0]


            if sum(normal**2) > 1.0e-32:
                normal = normal / numpy.sqrt(sum(normal**2))
            else:
                print normal
                print vec1, vec2
                raise Collision.CollisionException

            normal = normal * numpy.sign(numpy.dot(normal, (p-pa)))

            p = x + dt * (k - (1.0 + self.e) * normal * (numpy.dot(normal, k)))

            theta = abs(numpy.arcsin(numpy.dot(normal, (x-pa))
                                     / numpy.sqrt(numpy.dot(x - pa, x - pa))))

            coldat = []

            if any(v):
                vs = v + s * dt * f
            else:
                vs = self.v + s * dt * f

#            print 'Before', p0, p, vs

            par_col = copy.copy(self)
            par_col.p = x
            par_col.v = vs

            coldat.append(Collision.CollisionInfo(par_col, cell_index, theta,
                                                  self.t + s * dt, normal))
            vs += -(1.0 + self.e)* normal * numpy.dot(normal, vs)

#            print 'After V1:', pa, p, n, vs, f

            px, col, vo = self.collide(vs, (1 - s) * dt,
                                       v=vs, f=f, pa=x + 1.0e-10 * vs,
                                       level=level + 1)
            p = px + x + 1.0e-10 * vs

#            print 'After V2:', pa, p, n, vs, f

            if col:
                coldat += col

            return p - pa, coldat, vo

        return p - pa, None, v + dt * f - self.v

class ParticleBucket(object):
    """Class for a container for multiple Lagrangian particles."""

    def __init__(self, X, V, t=0, dt=1.0e-3, filename=None,
                 base_name='', U=None, GP=None, rho=2.5e3, g=numpy.zeros(3),
                 omega=numpy.zeros(3), diameter=40.e-6, boundary=None, e=0.99, tc=None):
        """Initialize the bucket 

        Args:
            X (float): Initial particle positions.
            V (float): Initial velocities
        """


        if tc:
            self.tc = tc
        else:
            self.tc = TemporalCache.TemporalCache(base_name)
        self.particles = []

        if U is None:
            U = [None for _ in range(X.shape[0])]
        if GP is None:
            GP = [None for _ in range(X.shape[0])]

        for pos, vel, fluid_vel, grad_p in zip(X, V, U, GP):
            self.particles.append(Particle(pos, vel, t, dt, tc=self.tc, u=fluid_vel,
                                           gp=grad_p, rho=rho, g=g, omega=omega,
                                           diameter=diameter, boundary=boundary, e=e))
        self.time = t
        self.pos = X
        self.vel = V
        self.fluid_vel = U
        self.grad_p = GP
        self.delta_t = dt
        self.boundary = boundary
        if filename:
            self.outfile = open(filename, 'w')

    def update(self):
        """ Update all the particles in the bucket to the next time level."""
        self.tc.range(self.time, self.time + self.delta_t)
        for part in self.particles:
            part.update()

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

    def run(self, time):
        """Drive particles forward until a given time."""
        while self.time < time:
            self.update()
            self.write()
        self.outfile.flush()
