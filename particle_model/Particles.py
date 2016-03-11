""" Baseline module for the package. Contains the main classes, particle and particle_bucket. """

from particle_model import IO
from particle_model import DragModels
from particle_model import Collision
from particle_model import System
from particle_model import Options
from particle_model import ParticleBase

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

class Particle(ParticleBase.ParticleBase):
    """Class representing a single Lagrangian particle with mass"""

    def __init__(self, data,
                 parameters=ParticleBase.PhysicalParticle(),
                 system=System.System(), **kwargs):

        super(Particle, self).__init__(*data, **kwargs)

        self.collisions = []
        self.parameters = parameters
        self.system = system
        self.solid_pressure_gradient=numpy.zeros(3)
        self.volume = self.parameters.get_volume()

    def update(self, delta_t=None, method=None):
        """ Update the state of the particle to the next time level."""
        if delta_t is not None:
            self.delta_t = delta_t
        if method=="AdamsBashforth2":
            self.update_ab2()
        else:
            self.update_rk4()

    def update_ab2(self):
        """Update the state of the particle to the next time level

        The method uses the Adams Bashforth second order method"""

        if self._old:

            kap = (self.vel, self.force(self.pos,
                                     self.vel,
                                     self.time), self.time)


            beta= 0.5*self.delta_t/(self.time-self._old[2])
            
            step, col, vel, col_vel = self.collide((1.0+beta)*self.vel-beta*self._old[0],
                                          self.delta_t,
                                          self.vel,
                                          force=(1.0+beta)*kap[1]-beta*self._old[1])

        else:
            ## reduced to using the Euler method for the first timestep:

            kap = (self.vel, self.force(self.pos,
                                     self.vel,
                                     self.time), self.time)

            step, col, vel, col_vel = self.collide(self.vel,
                                          self.delta_t,
                                          self.vel,
                                          force=kap[1])

        self.pos += step
        self.vel += vel
        if col:
            self.collisions += col

            kap = (self.vel, self.force(col[-1].pos+1.0e-8*col_vel,
                                     col_vel,
                                     col[-1].time+1.0e-8), col[-1].time+1.0e-8)

        self._old=kap
            
        self.time += self.delta_t

    def update_rk4(self):
        """Update the state of the particle to the next time level

        The method uses relatively simple RK4 time integration."""

        kap1 = (self.vel, self.force(self.pos,
                                     self.vel,
                                     self.time))

        step, col, vel, col_vel = self.collide(kap1[0], 0.5 * self.delta_t,
                                      self.vel,
                                      force=kap1[1])

        kap2 = (self.vel + 0.5*self.delta_t * kap1[1],
                self.force(self.pos + step,
                           self.vel + vel,
                           self.time + 0.5 * self.delta_t))

        step, col, vel, col_vel = self.collide(kap2[0], 0.5 * self.delta_t,
                                      vel=self.vel, force=kap2[1])
        kap3 = (self.vel+0.5 * self.delta_t * kap2[1],
                self.force(self.pos + step,
                           self.vel + vel,
                           self.time + 0.5*self.delta_t))

        step, col, vel, col_vel = self.collide(kap3[0], self.delta_t,
                                      vel=self.vel,
                                      force=kap3[1])
        kap4 = (self.vel + self.delta_t * kap3[1],
                self.force(self.pos + step,
                           self.vel + vel,
                           self.time + self.delta_t))


        step, col, vel, col_vel = self.collide((kap1[0] + 2.0 * (kap2[0] + kap3[0]) + kap4[0]) / 6.0,
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

        return (-grad_p / self.parameters.rho
                + self.parameters.drag(fluid_velocity,
                                       particle_velocity,
                                       diameter=self.parameters.diameter,
                                       rho=self.parameters.rho,
                                       rho_f=self.system.rho,
                                       fluid_viscosity=self.system.viscosity)
                / self.parameters.rho
                + self.coriolis_force(particle_velocity)
                + self.system.gravity
                + self.centrifugal_force(position)
                + self.solid_pressure_gradient / self.parameters.rho)

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
#        locator.FindCell(point)

        if cell_index==-1:
            cell_index = None

        return cell_index

    def picker(self, pos, time):
        """ Extract fluid velocity and pressure from .vtu files at correct time level"""

        def fpick(infile, locator,names):
            """ Extract fluid velocity and pressure from single .vtu file"""

            locator.BuildLocatorIfNeeded()

            cell_index = self.find_cell(locator, pos)
            cell = IO.get_linear_block(infile).GetCell(cell_index)
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

            data_u = IO.get_vector(infile, names[0], cell_index)
            data_p = IO.get_scalar(infile, names[1], cell_index)


            shape_funs = numpy.zeros(cell.GetNumberOfPoints())
            deriv_funs = numpy.zeros(dim*cell.GetNumberOfPoints())
            cell.InterpolateFunctions(upos[:3], shape_funs)
            cell.InterpolateDerivs(upos[:3], deriv_funs)

            rhs = data_p[1:dim+1]-data_p[0]

            mat = numpy.zeros((dim, dim))

            pts = numpy.array([cell.GetPoints().GetPoint(i) for i in range(dim+1)])
            for i in range(dim):
                mat[i, :] = pts[i+1, :dim] - pts[0, :dim]

#            mat=la.inv(mat)
            mat = invert(mat)

            out = numpy.dot(shape_funs, data_u)

            grad_p = numpy.zeros(3)
            grad_p[:dim] = numpy.dot(mat, rhs)

            return out, grad_p

        data, alpha, names = self.system.temporal_cache(time)

        vel0, grad_p0 = fpick(data[0][2], data[0][3],names[0])
        vel1, grad_p1 = fpick(data[1][2], data[1][3],names[1])

        vel0=numpy.append(vel0,numpy.zeros(3-vel0.size))
        vel1=numpy.append(vel1,numpy.zeros(3-vel1.size))


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
            return k * delta_t, None, vel - self.vel, None

        pos = pa+delta_t*k

        s = vtk.mutable(-1.0)
        x = [0.0, 0.0, 0.0]
        cell_index = vtk.mutable(0)

        bndl = self.system.boundary.bndl

        intersect = bndl.IntersectWithLine(pa, pos,
                                           1.0e-8, s,
                                           x, ARGV, ARGI, cell_index)

        if intersect:
            if self.system.boundary.bnd.GetCellData().HasArray('SurfaceIds'):
                if self.system.boundary.bnd.GetCellData().GetScalars('SurfaceIds').GetValue(cell_index) in self.system.boundary.outlet_ids:
                    return pos - pa, None, vel + delta_t * force - self.vel, None

            data, _, names = self.system.temporal_cache(self.time)
#            assert IO.test_in_cell(IO.get_linear_block(data[0][2]).GetCell(self.find_cell(data[0][3], pa)), pa) or sum((pa-x)**2)<1.0-10

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

            px, col, velo, dummy_vel = self.collide(vels, (1 - s) * delta_t,
                                         vel=vels, force=force,
                                         pa=x + 1.0e-9 * vels,
                                         level=level + 1)
            pos = px + x + 1.0e-9 * vels

            if col:
                coldat += col

            return pos - pa, coldat, velo, vels

        return pos - pa, None, vel + delta_t * force - self.vel, None

class ParticleBucket(object):
    """Class for a container for multiple Lagrangian particles."""

    def __init__(self, X, V, time=0, delta_t=1.0e-3, filename=None,
                 parameters=ParticleBase.PhysicalParticle(),
                 system=System.System()):
        """Initialize the bucket

        Args:
            X (float): Initial particle positions.
            V (float): Initial velocities
        """

        self.system = system

        ### pick only points which are actually in our test box
        live=system.in_system(X, time)
        X = X.compress(live, axis=0)
        V = V.compress(live, axis=0)

        self.fluid_velocity = numpy.zeros(X.shape)
        self.grad_p = numpy.zeros(X.shape)

        self.particles = []
        self.dead_particles = []
        self.parameters = parameters
        for _, (dummy_pos, dummy_vel) in enumerate(zip(X, V)):
            self.particles.append(Particle((dummy_pos, dummy_vel, time, delta_t),
                                           system=self.system,
                                           parameters=parameters.randomize()))
            if self.system.temporal_cache:
                dummy_u, dummy_p = self.particles[-1].get_fluid_properties()
                self.fluid_velocity[_, :len(dummy_u)]=dummy_u
                self.grad_p[_, :] = dummy_p
        self.time = time
        self.delta_t = delta_t
        self.pos = X
        self.vel = V
        self.solid_pressure_gradient = numpy.zeros((len(self.particles),3))
        for particle, gsp in zip(self.particles,self.solid_pressure_gradient):
            particle.solid_pressure_gradient=gsp
            
        if filename:
            self.outfile = open(filename, 'w')

    def update(self, delta_t=None,*args,**kwargs):
        if delta_t is not None:
            self.delta_t = delta_t
        """ Update all the particles in the bucket to the next time level."""
        self.system.temporal_cache.range(self.time, self.time + self.delta_t)
        live=self.system.in_system(self.pos, self.time)
        for k, part in enumerate(self.particles):
            if live[k]:
                part.update(self.delta_t,*args,**kwargs)
            else:
                self.dead_particles.append(part)
                self.particles.remove(part)
        self.insert_particles()
        if self.system.temporal_cache:
            self.fluid_velocity=numpy.empty((len(self.particles),3),float)
            self.grad_p=numpy.empty((len(self.particles),3),float)
            self.pos=numpy.empty((len(self.particles),3),float)
            self.vel=numpy.empty((len(self.particles),3),float)
            for k, part in enumerate(self.particles):
                self.pos[k,:]=part.pos
                self.vel[k,:]=part.vel
                self.fluid_velocity[k, :], self.grad_p[k, :] = \
                        part.get_fluid_properties()

        self.time += self.delta_t

    def insert_particles(self):
        """Deal with particle insertion"""

        for inlet in self.system.boundary.inlets:
            n_par = inlet.get_number_of_insertions(self.time,
                                                   self.delta_t)
            if n_par == 0:
                continue
            weights = inlet.cum_weight(self.time+0.5*self.delta_t,
                                       self.system.boundary.bnd)
            for i in range(n_par):
                prob = numpy.random.random()
                time =  self.time+prob*self.delta_t 
                pos =inlet.select_point(time, weights,
                                        self.system.boundary.bnd)
                vel = numpy.array(inlet.velocity(pos, time))

                par = Particle((pos, vel, time,
                                (1.0-prob)*self.delta_t),
                               system=self.system,
                               parameters=self.parameters.randomize())

                par.update()
                par.time = self.time + self.delta_t
                par.delta_t = self.delta_t

                self.particles.append(par)


    def collisions(self):
        """Collect all collisions felt by particles in the bucket"""
        return [i for i in itertools.chain(*[p.collisions for p in self.particles+self.dead_particles])]


    def write(self):
        """Write timelevel data to file."""

        self.outfile.write('%f'%self.time)

        for pos in self.pos.ravel():
            self.outfile.write(' %f'%pos)

        for vel in self.vel.ravel():
            self.outfile.write(' %f'%vel)

        self.outfile.write('\n')

    def set_solid_pressure_gradient(self,solid_pressure_gradient):
        self.solid_pressure_gradient[:,:]=solid_pressure_gradient
        for particle, gsp in zip(self.particles,self.solid_pressure_gradient):
            particle.solid_pressure_gradient=gsp

    def run(self, time, delta_t=None, write=True,*args,**kwargs):
        """Drive particles forward until a given time."""
        while self.time - time < -1.0e-15:
            self.update(delta_t,*args,**kwargs)
            if write:
                self.write()
        if write:
            self.outfile.flush()
