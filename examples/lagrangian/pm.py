"""Example of passive Lagrangian particles in a 2d double gyre hitting wall."""

import sys

import particle_model as pm
import numpy

# problem parameters
if len(sys.argv)<2:
    S='40'
else:
    S = sys.argv[1]

N = 2000  # number of particles

X = 0.5+0.25*(numpy.random.random((N, 3))-0.5) #initial positions
X[:, 0] += 0.2
X[:, 2] = 0

V = numpy.zeros((N, 3)) # initial velocities, not used for passive Lagrangian

NAME = 'gyre%sx%s'%(S,S) # base name for .vtus

TEMP_CACHE = pm.TemporalCache.TemporalCache(NAME, 
                                            velocity_name="Velocity", 
                                            pressure_name=None,
                                            time_name="Time")

MESH = pm.IO.GmshMesh()
MESH.read('Structured%sx%s.msh'%(S,S))
BOUNDARY_MESH = pm.IO.make_boundary_from_msh(MESH)
INLET = pm.Options.Inlet(surface_ids=[1], insertion_rate=1.0e5,
                         velocity=lambda X,t,: (-2*numpy.pi*numpy.cos(2*numpy.pi*X[1]),0.0,0.0),
                         pdf= lambda X,t : max(0,-numpy.cos(2*numpy.pi*X[1]))) # setup inlet, inlet velocity and pdf
BOUNDARY = pm.IO.BoundaryData(BOUNDARY_MESH, inlets=[INLET])
SYSTEM = pm.System.System(BOUNDARY, temporal_cache=TEMP_CACHE, outlet_ids=[1], inlets=[INLET]) #outlets allow particles to leave

PAR = pm.ParticleBase.PhysicalParticle(diameter=0.0)

PB = pm.Particles.ParticleBucket(X, V, time=0.0, delta_t=1.0e-3,
                                 system=SYSTEM,
                                 parameters=PAR) #time is start time, delta_t is timestep

for i, level in enumerate(TEMP_CACHE):

    print 'time', level[0]

    PB.run(level[0], write=False,method="AdamsBashforth2")
    pm.IO.write_level_to_polydata(PB, i, NAME)

    print 'min, max: pos_x', PB.pos_as_array()[:, 0].min(), PB.pos_as_array()[:, 0].max()
    print 'min, max: pos_y', PB.pos_as_array()[:, 1].min(), PB.pos_as_array()[:, 1].max()
    print 'min, max: vel_x', PB.vel_as_array()[:, 0].min(), PB.vel_as_array()[:, 0].max()
    print 'min, max: vel_y', PB.vel_as_array()[:, 1].min(), PB.vel_as_array()[:, 1].max()
