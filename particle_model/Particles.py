""" Baseline module for the package. Contains the main classes, particle and particle_bucket. """

# standard imports
import itertools
import copy

import numpy
import vtk

from particle_model.Debug import profile, logger
from particle_model import Math
from particle_model import IO
from particle_model import Collision
from particle_model import System
from particle_model import ParticleBase
from particle_model import Parallel
from particle_model import TemporalCache
from particle_model import Timestepping
from particle_model import vtk_extras

LEVEL = 0
ZERO = numpy.zeros(3)

class Particle(ParticleBase.ParticleBase):
    """Class representing a single Lagrangian particle with mass"""

    def __init__(self, data,
                 parameters=ParticleBase.PhysicalParticle(),
                 system=System.System(), **kwargs):

        super(Particle, self).__init__(*data, **kwargs)

        self.collisions = []
        self.parameters = parameters
        self.pure_lagrangian = self.parameters.pure_lagrangian()
        self.system = system
        self.solid_pressure_gradient = numpy.zeros(3)
        self.volume = self.parameters.get_volume()

        self.pos_callbacks = kwargs.get('pos_callbacks', [])
        self.vel_callbacks = kwargs.get('vel_callbacks', [])

    def __repr__(self):
        return "Particle((%r, %r, %r, %r, %r) , %r, %r)"%(self.pos,
                                                          self.vel,
                                                          self.time,
                                                          self.delta_t,
                                                          self._hash,
                                                          self.parameters,
                                                          self.system)

    def copy(self):
        """ Create a (mixed) copy of the particle."""
        par = Particle((self.pos, self.vel, self.time, self.delta_t),
                       parameters=self.parameters, system=None)
        par.set_hash(self._hash)
        par._old = copy.deepcopy(self._old)
        return par

    def update(self, delta_t=None, method="AdamsBashforth2"):
        """ Update the state of the particle to the next time level."""
        if delta_t is not None:
            self.delta_t = delta_t
        try:
            Timestepping.methods[method](self)
        except KeyError:
            logger.warning("Timestepping method %s unknown, using AdamsBashforth2.",
                           method)
            Timestepping.methods["AdamsBashforth2"](self)
        except Collision.CollisionException as col:
            # wall collision occurred.
            self.vel = Collision.rebound_velocity(self, col.vel, numpy.zeros(3), col.info.normal, col.info.cell)
            self.pos = col.pos+1.0e-10*col.info.normal
            self.time += col.delta_t
            col.info.time = self.time
            self.collisions += [col.info]
            self._old = []
            self.update(self.delta_t-col.delta_t, method)
        except Collision.OutletException as col:
            self.exited = True


    def drag_coefficient(self, position, particle_velocity, time, nearest=False):
        """ Get particle drag coefficent for specified position and velocity. """

        fluid_velocity = self.picker(position, time, nearest=nearest)[0]

        if fluid_velocity is not None:
            drag = self.parameters.drag_coefficient(fluid_velocity,
                                                    particle_velocity,
                                                    diameter=self.parameters.diameter,
                                                    rho=self.parameters.rho,
                                                    rho_f=self.system.rho,
                                                    fluid_viscosity=self.system.viscosity)
        else:
            drag = 0.0

        return drag/self.parameters.rho, fluid_velocity

    def get_fluid_properties(self):
        """ Get the fluid velocity and pressure gradient at the particle
        location."""
        return self.picker(self.pos, self.time)

    @profile
    def force(self, position, particle_velocity, time, drag=True):
        """Calculate the sum of the forces on the particle.

        Args:
            p (float): Location at which forcing is evaluated.
            v (float): Particle velocity at which forcing is evaluated.
            t (float): Time at which particle is evaluated.
        """

#        if collision:
#            raise collisionException


        if self.pure_lagrangian:
            return ZERO

        fluid_velocity, grad_p = self.picker(position, time)

        if drag and (fluid_velocity is not None):
            drag_force = self.parameters.drag(fluid_velocity,
                                              particle_velocity,
                                              diameter=self.parameters.diameter,
                                              rho=self.parameters.rho,
                                              rho_f=self.system.rho,
                                              fluid_viscosity=self.system.viscosity)
        else:
            drag_force = 0.0
#            grad_p = numpy.zeros(3)
            fluid_velocity, grad_p = self.picker(position, time)

        if grad_p is None:
            grad_p = numpy.zeros(3)
            drag_force = 0.0

#        try:
        return (-1.0*grad_p / self.parameters.rho
                + drag_force/ self.parameters.rho
                + self.coriolis_force(particle_velocity)
                + self.system.gravity
                + self.centrifugal_force(position)
                + self.solid_pressure_gradient / self.parameters.rho)
#        except:
#            IPython.embed()

    def coriolis_force(self, particle_velocity):
        """ Return Coriolis force on particle."""
        return -2.0 * Math.cross(self.system.omega, particle_velocity)

    def centrifugal_force(self, position):
        """ Return centrifugal force on particle"""
        return - (numpy.dot(self.system.omega, position) * self.system.omega
                  - numpy.dot(self.system.omega, self.system.omega) * position)

    def find_cell(self, locator, point=None):
        """ Use vtk rountines to find cell/element containing the point."""
        if point is None:
            point = self.pos
        cell_index = vtk.mutable(0)

        cell_index, pcoords = vtk_extras.FindCell(locator, point)

        if cell_index == -1:
            cell_index = None

        return cell_index, pcoords

    @profile
    def _fpick(self, pos, infile, picker, names, nearest=False):
        """ Extract fluid velocity and pressure from single .vtu file"""

        if picker.cell_index is None and not nearest:
            return None, None

        # picker defaults to velocity
        if nearest:
            out = picker.nearest(pos)
        else:
            out = picker(pos)

        if names[1] and picker.cell_index:
            data_p = IO.get_scalar(infile,
                                   self.system.temporal_cache.get(infile, names[1]),
                                   names[1], picker.cell_index)
        else:
            data_p = None

        if data_p is not None:
            dim = picker.cell.GetCellDimension()
            pts = numpy.array([picker.cell.GetPoints().GetPoint(i) for i in range(picker.cell.GetNumberOfPoints())])
            if (pts.shape[0] == dim+1 
                or (pts.shape[0] == 6 and dim == 2)
                or (pts.shape[0] == 10 and dim == 3)):
                grad_p = Math.grad(data_p, pts, dim)
            else:
                grad_p = numpy.zeros(3)
                picker.cell.Derivatives(0, picker.pcoords, data_p, 1, grad_p)
                
        else:
            grad_p = ZERO

        if len(names) == 3:
            vdata = self.system.temporal_cache.get(infile, names[2])
            gout = picker(pos, vdata)
            return out, grad_p, gout
        #otherwise
        return out, grad_p


    @profile
    def picker(self, pos, time, gvel_out=None, nearest=False):
        """ Extract fluid velocity and pressure from .vtu files at correct time level"""

        data, alpha, names = self.system.temporal_cache(time)

        TemporalCache.PICKERS[0].name = names[0][0]
        TemporalCache.PICKERS[0].grid = IO.get_block(data[0][2], names[0][0])
        TemporalCache.PICKERS[0].locator = data[0][3]
        TemporalCache.PICKERS[0].pos = numpy.asarray(pos, float)
        TemporalCache.PICKERS[1].name = names[1][0]
        TemporalCache.PICKERS[1].grid = IO.get_block(data[1][2], names[1][0])
        TemporalCache.PICKERS[1].locator = data[1][3]
        TemporalCache.PICKERS[1].pos = numpy.asarray(pos, float)

        if len(names[0]) == 3:
            vel0, grad_p0, gvel_0 = self._fpick(pos, data[0][2],
                                                TemporalCache.PICKERS[0], names[0],
                                                nearest)
            vel1, grad_p1, gvel_1 = self._fpick(pos, data[1][2],
                                                TemporalCache.PICKERS[1], names[1],
                                                nearest)
        else:
            vel0, grad_p0 = self._fpick(pos, data[0][2],
                                        TemporalCache.PICKERS[0], names[0],
                                        nearest)
            vel1, grad_p1 = self._fpick(pos, data[1][2],
                                        TemporalCache.PICKERS[1], names[1],
                                        nearest)

        if vel0 is None or vel1 is None:
            return None, None

        if gvel_out:
            gvel_out = (1.0-alpha)*gvel_0+alpha*gvel_1

        return ((1.0-alpha) * vel0 + alpha * vel1,
                (1.0 - alpha) * grad_p0 + alpha * grad_p1)

    def _check_remapping(self, pos_1, pos_0, vel_0, delta_t):
        """Test for periodic/remapped boundaries"""

        ### this finds the point of boundary intersection
        ### pos_i = pos_0 + t_val*(pos_1-pos_0)

        intersect, pos_i, t_val, cell_index, pcoords = self.system.boundary.test_intersection(pos_0, pos_1)

        if intersect and cell_index >= 0:
            surface_id = self.system.boundary.get_surface_id(cell_index)
            if surface_id is not None:
                if surface_id in self.system.boundary.mapped_ids:
                    pos_o, vel_o = self.system.boundary.mapped_ids[surface_id](pos_i, vel_0)

                    pos_f = pos_o+(1.0-t_val)*delta_t*vel_o

                    vel_i, _ = self.picker(pos_f, self.time+delta_t)

                    par_col = copy.copy(self)
                    if self.system.boundary.dist:
                        cell = self.system.boundary.bnd.GetCell(cell_index)
                        par_col.pos = IO.get_real_x(cell, pcoords)
                    else:
                        par_col.pos = pos_i
                    par_col.vel = vel_0
                    par_col.time = self.time + t_val * delta_t

                    return pos_f, vel_i
                elif surface_id in self.system.boundary.outlet_ids: 
                    raise Collision.OutletException(pos_1, vel_0)
                else:
                    # This is a "reflecting" boundary.
                    angle, normal = Collision.collision_angle(self, pos_0, pos_i,
                                        cell_index)
                    return pos_1 - 2.0*numpy.dot(pos_1-pos_i, normal)*normal, vel_0
                    

        #otherwise
        
        return pos_0, vel_0

    def check_collision_full(self, pos_1, pos_0, vel_1, vel_0, delta_t, drag):
        """ Check for particle-wall collision.

        Args:
            pos_1, pos_1  (float): Initial and final displacement
            vel_1, vel_0 (float): Initial and final velocity
            delta_t (float): timestep
            drag (boolean): do drag implicitly.

        Return final position and velocity.
        """

        if self.pure_lagrangian:
            fvel = self.picker(pos_1, self.time+delta_t)[0]
            if fvel is None:
                return self._check_remapping(pos_1, pos_0, vel_0, delta_t)

            for cback in self.vel_callbacks:
                fvel += delta_t*cback(pos_1, fvel, self.time+delta_t, delta_t)

            return pos_1, fvel

        intersect, pos_i, t_val, cell_index, pcoords = self.system.boundary.test_intersection(pos_0, pos_1)

        if intersect and cell_index >= 0:
            surface_id = self.system.boundary.get_surface_id(cell_index)
            if surface_id is not None:
                if surface_id in self.system.boundary.outlet_ids:
                    raise Collision.OutletException(pos_1, vel_1)
                elif surface_id in self.system.boundary.mapped_ids:
                    raise Collision.MappedBoundaryException(self.system.boundary.mapped_ids[surface_id])
            #otherwise
            raise Collision.CollisionException(self, pos_i, cell_index,
                                               t_val*delta_t)

        # no collisions
        if drag:
            return pos_1, self._drag_update(pos_1, vel_0, vel_1,
                                            self.time+delta_t, delta_t)
        return pos_1, vel_1

    def _drag_update(self, pos, vel_0, vel_1, time, delta_t):
        """ Solve for drag term implictly."""

        c_d, fvel = self.drag_coefficient(pos, vel_1, time)
        if fvel is None:
            return vel_1/(1.0+delta_t*c_d)
        return (vel_1+delta_t*c_d*fvel)/(1.0+delta_t*c_d)

class ParticleBucket(object):
    """Class for a container for multiple Lagrangian particles."""

    def __init__(self, X, V, time=0, delta_t=1.0e-3,
                 parameters=ParticleBase.PhysicalParticle(),
                 system=System.System(),
                 field_data=None, online=True, **kwargs):
        """Initialize the bucket

        Args:
            X (float): Initial particle positions.
            V (float): Initial velocities
        """

        logger.info("Initializing ParticleBucket")

        field_data = field_data or {}

        self.system = system

        ### pick only points which are actually in our test box
        live = system.in_system(X, len(X), time)
        X = X.compress(live, axis=0)
        V = V.compress(live, axis=0)

        self.particles = []
        self.dead_particles = []
        self.stuck_particles = []
        self.parameters = parameters
        for _, (dummy_pos, dummy_vel) in enumerate(zip(X, V)):
            par = Particle((dummy_pos, dummy_vel, time, delta_t),
                           system=self.system,
                           parameters=parameters.randomize(),
                           **kwargs)
            if par.pure_lagrangian:
                par.vel = par.picker(par.pos, time)[0]
            for name, value in field_data.items():
                par.fields[name] = copy.deepcopy(value[_])
            if self.system.temporal_cache:
                dummy_u, dummy_p = par.get_fluid_properties()
                if dummy_u is not None:
                    self.particles.append(par)
            else:
                self.particles.append(par)
        self.time = time
        self.delta_t = delta_t
        self._online = online
        self.solid_pressure_gradient = numpy.zeros((len(self.particles), 3))
        for particle, gsp in zip(self, self.solid_pressure_gradient):
            particle.solid_pressure_gradient = gsp

        self.redistribute()

    def __len__(self):
        return len(self.particles)

    def __bool__(self):
        return bool(self.particles)

    def __contains__(self, par):
        return par in self.particles

    def __nonzero__(self):
        return self.__bool__()

    def __iter__(self):
        for particle in self.particles:
            yield particle

    def pos(self):
        """Generator function for particle positions."""
        for part in self:
            yield part.pos

    def pos_as_array(self):
        """Particle positions as numpy array."""
        out = numpy.empty((len(self), 3), float)
        for _, pos in enumerate(self.pos()):
            out[_, :] = pos[:]
        return out

    def vel(self):
        """Generator function for particle velocities."""
        for part in self:
            yield part.vel

    def vel_as_array(self):
        """Particle velocities as numpy array."""
        out = numpy.empty((len(self), 3), float)
        for _, vel in enumerate(self.vel()):
            out[_, :] = vel[:]
        return out

    @profile
    def update(self, delta_t=None, *args, **kwargs):
        """ Update all the particles in the bucket to the next time level."""

        logger.info("In ParticleBucket.Update: %d particles", len(self.particles))

        # redistribute particles to partitions in case of parallel adaptivity
        if Parallel.is_parallel():
            self.redistribute()
        # reset the particle timestep
        if delta_t is not None:
            self.delta_t = delta_t
        self.system.temporal_cache.range(self.time, self.time + self.delta_t)
        live = self.system.in_system(self.pos(), len(self), self.time)
        _ = []
        for k, part in enumerate(self):
            if live[k] and not hasattr(part, "exited"):
                try: 
                #if particle updates fails e.g. max recursion depth reached
                    part.update(self.delta_t, *args, **kwargs)
                except RuntimeError: 
                #remove as a stuck particle
                    self.stuck_particles.append(part)
                    _.append(part)
            else:
                self.dead_particles.append(part)
                _.append(part)
        for part in _:
            self.particles.remove(part)
        self.redistribute()
        self.insert_particles(*args, **kwargs)
        self.time += self.delta_t
        for part in self:
            part.time = self.time

    def redistribute(self):
        """ In parallel, redistrbute particles to their owner process."""
        if self._online and Parallel.is_parallel():
            logger.debug("%d particles before redistribution", len(self.particles))
            self.particles = Parallel.distribute_particles(self.particles,
                                                           self.system)

            logger.debug("%d particles after redistribution", len(self))

    def insert_particles(self, *args, **kwargs):
        """Deal with particle insertion"""

        for inlet in self.system.boundary.inlets:
            if Parallel.is_parallel():
                n_par = inlet.get_number_of_insertions(self.time,
                                                       self.delta_t)
#                if Parallel.get_rank() == 0:
#                    n_par_0 = inlet.get_number_of_insertions(self.time,
#                                                           self.delta_t)
#                    for i in range(n_par_0):
#                        prob = numpy.random.random()
            else:
                n_par = inlet.get_number_of_insertions(self.time,
                                                       self.delta_t)
            if n_par == 0:
                continue
            weights = inlet.cum_weight(self.time+0.5*self.delta_t,
                                       self.system.boundary.bnd,
                                       self.system.temporal_cache)
            if weights:
                for i in range(n_par):
                    prob = numpy.random.random()
                    time = self.time+prob*self.delta_t
                    pos = inlet.select_point(time, weights,
                                             self.system.boundary.bnd)
                    if inlet.velocity:
                        vel = numpy.array(inlet.velocity(pos, time))
                    else:
                        vel = numpy.zeros(3)
                        fvel = self.system.temporal_cache.get_velocity(pos, time)
                        vel[:len(fvel)] = fvel

                    ## update position by fractional timestep
                    pos = pos + vel*(1-prob)*self.delta_t

                    data, alpha, names = self.system.temporal_cache(time)
                    cell_id, pcoords = vtk_extras.FindCell(data[0][3], pos)

                    if cell_id == -1:
                        continue

                    par = Particle((pos, vel, time,
                                    (1.0-prob)*self.delta_t),
                                   system=self.system,
                                   parameters=self.parameters.randomize(),
                                   **inlet.kwargs)

                    par.delta_t = self.delta_t

                    par.fields["InsertionTime"] = time
                    self.particles.append(par)

    def collisions(self):
        """Collect all collisions felt by particles in the bucket"""
        return [i for i in itertools.chain(*[p.collisions for p in self.particles+self.dead_particles])]

    def set_solid_pressure_gradient(self, solid_pressure_gradient):
        self.solid_pressure_gradient = solid_pressure_gradient
        for particle, gsp in zip(self.particles, self.solid_pressure_gradient):
            particle.solid_pressure_gradient = gsp

    def run(self, time, delta_t=None, write=False, *args, **kwargs):
        """Drive particles forward until a given time."""
        while time-self.time > 1.0e-6*(delta_t or self.delta_t):
            self.update(delta_t, *args, **kwargs)
            if write:
                global LEVEL
                IO.write_level_to_polydata(bucket=self, level=LEVEL, basename="dump.vtp")
                LEVEL += 1
