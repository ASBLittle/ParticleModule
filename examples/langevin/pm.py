"""Example of schochastic dispersion of passive Lagrangian particles in a 2d quadruple gyre."""

import sys
import numpy
import numpy.random

## Uncomment this to profile code
#from fluidity_tools import Profile
#try:
#    import builtins
#except ImportError:
#    import __builtin__ as builtins
#builtins.__dict__['profile'] = Profile.make_line_profiler()

import particle_model as pm

# problem parameters

kap = float(sys.argv[1])
if len(sys.argv) < 3:
    S = '40'
else:
    S = sys.argv[2]

def h(pos):
    """ a depth term"""
    return 1.0-((pos[0]-0.5)**2+(pos[1]-0.5)**2)

def grad_h(pos):
    """ a depth gradient term"""
    return numpy.array((-2.0*(pos[0]-0.5), -2.0*(pos[1]-0.5), 0.0))

def kappa(pos):
    """ a diffusion coeffient """
    return kap*(1.0-((pos[0]-0.5)**2+(pos[1]-0.5)**2))

def grad_kappa(pos):
    """ a diffusion coeffient gradient """
    return kap*numpy.array((-2.0*(pos[0]-0.5), -2.0*(pos[1]-0.5), 0.0))

def conservation_correction(pos, vel, time, delta_t):
    """ q.v. Spangol et al."""
    return (kappa(pos)/h(pos)*grad_h(pos)+grad_kappa(pos))

rwalk = pm.Stochastic.nonconstant_gaussian_drift_term(kappa=kappa, dim=2)
# rwalk is a function with signature func(pos, vel, time, delta_t)
    



N = 500  # number of particles

X = 0.5+0.25*(numpy.random.random((N, 3))-0.5) #initial positions
X[:, 2] = 0

V = numpy.zeros((N, 3)) # initial velocities, not used for passive Lagrangian

NAME = 'gyre%sx%s'%(S, S) # base name for .vtus

TEMP_CACHE = pm.TemporalCache.TemporalCache(NAME,#   use NAME+'.pvd' for .pvd reader 
                                            velocity_name="Velocity", 
                                            pressure_name=None,
                                            time_name="Time",
                                            online=False)

for k, x in enumerate(X):
    V[k, :] = TEMP_CACHE.get_velocity(x, 0.0)


MESH = pm.IO.GmshMesh()
MESH.read('Structured%sx%s.msh'%(S, S))
BOUNDARY_MESH = pm.IO.make_boundary_from_msh(MESH)
BOUNDARY = pm.IO.BoundaryData(BOUNDARY_MESH, inlets=[])
SYSTEM = pm.System.System(BOUNDARY, temporal_cache=TEMP_CACHE, outlet_ids=[], inlets=[]) #outlets allow particles to leave

PAR = pm.ParticleBase.PhysicalParticle(diameter=0.0)


orig_gyre = numpy.zeros((N, 1))

for _ in range(N):
    orig_gyre[_, 0] = 2.0*(X[_, 0]>0.5)+(X[_, 1]>0.5)

FIELDS = {"InsertionTime": numpy.zeros((N, 1)),
          "OriginalGyre": orig_gyre }

# buckets hold collections of particles
# drive the system from there

PB = pm.Particles.ParticleBucket(X, V, time=0.0, delta_t=5.0e-4,
                                 system=SYSTEM,
                                 parameters=PAR,
                                 field_data=FIELDS,
                                 online=False,
                                 pos_callbacks=[rwalk,conservation_correction])
#time is start time, delta_t is timestep. Setting online only matters in parallel

#example of setting field data in script, iterate over bucket

for particle in PB:
    particle.fields["ExampleLevel"] = 1.0

PD = pm.IO.PolyData(NAME+'_trajectories', 
                    {"InsertionTime":1,
                     "OriginalGyre":1,
                     "ExampleLevel":1}) # This holds trajectory information
# output format is dictionary, key is name, value is length

PD.append_data(PB) # Store initial particle positions

for i, cache in enumerate(TEMP_CACHE):

    print('time', cache[0])

    # call which updates the particles
    PB.run(time=cache[0], write=False, method="AdamsBashforth2")

    # lets explicitly insert a new particle
    xx = numpy.array((0.3, 0.1, 0.0)) # Space is always 3d: in 2d z component is zero
    vv = TEMP_CACHE.get_velocity(xx, PB.time) # Same for velocity.
    part = pm.Particles.Particle((xx, vv, PB.time, PB.delta_t),
                                 system=SYSTEM,
                                 parameters=PAR,
                                 pos_callbacks=[rwalk, conservation_correction]) # Make a new particle
    part.fields["InsertionTime"] = part.time # set its insertion time
    part.fields["OriginalGyre"] = 0.0
    if pm.Parallel.get_rank() == 0:
        PB.particles.append(part) # Stick it in the bucket
    
    #example of setting field data inside script
    for particle in PB:
        t = particle.time
        t_0 = particle.fields["InsertionTime"]
        particle.fields["ExampleLevel"] = numpy.exp(-(t-t_0)**2)

    pm.IO.write_level_to_csv(bucket=PB, level=i, basename=NAME+'_%s'%kap,
                                  field_data={"InsertionTime":1,
                                              "OriginalGyre":1,
                                              "ExampleLevel":1}) # Dump just this timelevel
    pm.IO.write_level_to_polydata(bucket=PB, level=i, basename=NAME+'_%s'%kap,
                                  field_data={"InsertionTime":1,
                                              "OriginalGyre":1,
                                              "ExampleLevel":1}) # Dump just this timelevel
    PD.append_data(PB)

    print('min, max: pos_x', PB.pos_as_array()[:, 0].min(), PB.pos_as_array()[:, 0].max())
    print('min, max: pos_y', PB.pos_as_array()[:, 1].min(), PB.pos_as_array()[:, 1].max())
    print('min, max: vel_x', PB.vel_as_array()[:, 0].min(), PB.vel_as_array()[:, 0].max())
    print('min, max: vel_y', PB.vel_as_array()[:, 1].min(), PB.vel_as_array()[:, 1].max())

PD.write() # write trajectories
