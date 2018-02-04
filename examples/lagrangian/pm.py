"""Example of passive Lagrangian particles in a 2d double gyre hitting wall."""

import sys

import particle_model as pm
import numpy

# problem parameters
if len(sys.argv)<2:
    S='40'
else:
    S = sys.argv[1]

N = 500  # number of particles

X = 0.5+0.25*(numpy.random.random((N, 3))-0.5) #initial positions
X[:, 0] += 0.2
X[:, 2] = 0

V = numpy.zeros((N, 3)) # initial velocities, not used for passive Lagrangian

NAME = 'gyre%sx%s'%(S,S) # base name for .vtus

TEMP_CACHE = pm.TemporalCache.TemporalCache(NAME, 
                                            velocity_name="Velocity", 
                                            pressure_name=None,
                                            time_name="Time",
                                            online=False)

MESH = pm.IO.GmshMesh()
MESH.read('Structured%sx%s.msh'%(S,S))
BOUNDARY_MESH = pm.IO.make_boundary_from_msh(MESH)
INLET = pm.Options.Inlet(surface_ids=[1], insertion_rate=5.0e3,
                         velocity=lambda X,t,: (-2*numpy.pi*numpy.cos(2*numpy.pi*X[1]),0.0,0.0),
                         pdf= lambda X,t : max(0,-numpy.cos(2*numpy.pi*X[1]))) # setup inlet, inlet velocity and pdf
BOUNDARY = pm.IO.BoundaryData(BOUNDARY_MESH, inlets=[INLET])
SYSTEM = pm.System.System(BOUNDARY, temporal_cache=TEMP_CACHE, outlet_ids=[1], inlets=[INLET]) #outlets allow particles to leave

PAR = pm.ParticleBase.PhysicalParticle(diameter=0.0)

FIELDS = { "InsertionTime": numpy.zeros((N,1)) }

# buckets hold collections of particles
# drive the system from there

PB = pm.Particles.ParticleBucket(X, V, time=0.0, delta_t=5.0e-3,
                                 system=SYSTEM,
                                 parameters=PAR,
                                 field_data=FIELDS,
                                 online=False) 
#time is start time, delta_t is timestep. Setting online only matters in parallel

#example of setting field data in script, iterate over bucket

for particle in PB:
    particle.fields["ExampleLevel"] = 1.0

PD = pm.IO.PolyData(NAME+'_trajectories', 
                    { "InsertionTime":1,
                      "ExampleLevel":1}) # This holds trajectory information
# output format is dictionary, key is name, value is length
 
PD.append_data(PB) # Store initial particle positions

for i, cache in enumerate(TEMP_CACHE):

    print 'time', cache[0]

    # call which updates the particles
    PB.run(time=cache[0], write=False, method="AdamsBashforth2")

    # lets explicitly insert a new particle
    xx = numpy.array((0.3, 0.1, 0.0)) # Space is always 3d: in 2d z component is zero
    vv = numpy.array((0.0, 0.0, 0.0)) # Same for velocity.
    part = pm.Particles.Particle((xx, vv, PB.time, PB.delta_t),
                       system=SYSTEM,
                       parameters=PAR) # Make a new particle
    part.fields["InsertionTime"] = part.time # set its insertion time
    if pm.Parallel.get_rank()==0: PB.particles.append(part) # Stick it in the bucket
    
    #example of setting field data inside script
    for particle in PB:
        t = particle.time
        t_0 = particle.fields["InsertionTime"]
        particle.fields["ExampleLevel"] = numpy.exp(-(t-t_0)**2)

    pm.IO.write_level_to_polydata(bucket=PB, level=i, basename=NAME,
                                  field_data={ "InsertionTime":1,
                                               "ExampleLevel":1}) # Dump just this timelevel
    PD.append_data(PB)

    print 'min, max: pos_x', PB.pos_as_array()[:, 0].min(), PB.pos_as_array()[:, 0].max()
    print 'min, max: pos_y', PB.pos_as_array()[:, 1].min(), PB.pos_as_array()[:, 1].max()
    print 'min, max: vel_x', PB.vel_as_array()[:, 0].min(), PB.vel_as_array()[:, 0].max()
    print 'min, max: vel_y', PB.vel_as_array()[:, 1].min(), PB.vel_as_array()[:, 1].max()

PD.write() # write trajectories
