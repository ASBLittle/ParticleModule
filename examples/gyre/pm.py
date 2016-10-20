"""Example of particles in a 2d double gyre hitting wall."""
import particle_model as pm
import numpy

S = '40'

N = 2000

X = 0.5+0.25*(numpy.random.random((N, 3))-0.5)
X[:, 0] += 0.2
X[:, 2] = 0

V = numpy.zeros((N, 3))
V[:, 0] = -2.0*numpy.pi*numpy.sin(numpy.pi*X[:, 0])*numpy.cos(2.0*numpy.pi*X[:, 1])
V[:, 1] = numpy.pi*numpy.cos(numpy.pi*X[:, 0])*numpy.sin(2.0*numpy.pi*X[:, 1])

NAME = 'gyre'

TEMP_CACHE = pm.TemporalCache.TemporalCache(NAME)
BOUNDARY = pm.IO.BoundaryData('Gyre_boundary.vtu')
SYSTEM = pm.System.System(BOUNDARY, temporal_cache=TEMP_CACHE)
PAR = pm.ParticleBase.PhysicalParticle(diameter=1e-3)

PB = pm.Particles.ParticleBucket(X, V, 0.0, 1.0e-3,
                                 system=SYSTEM,
                                 parameters=PAR)

TEMP_CACHE.data[1][0] = 100.0

PD = pm.IO.PolyData(NAME+'.vtp')
PD.append_data(PB)
gsp = pm.IO.write_level_to_polydata(PB, 0, NAME, do_average=True)
PB.set_solid_pressure_gradient(gsp)

for i in range(300):
    print PB.time
    print 'min, max: pos_x', PB.pos[:, 0].ravel().min(), PB.pos[:, 0].ravel().max()
    print 'min, max: vel_x', PB.vel[:, 0].ravel().min(), PB.vel[:, 0].ravel().max()
    PB.update()
    PD.append_data(PB)
    gsp = pm.IO.write_level_to_polydata(PB, i+1, NAME, do_average=True)
    PB.set_solid_pressure_gradient(gsp)
#    ugrid = pm.IO.point_average(TEMP_CACHE.data[0][2], PB)
#    pm.IO.write_to_file(ugrid, "%s_out_%d.vtu"%('my_gyre', i))
#    pm.IO.write_level_to_ugrid(PB, i+1, NAME,
#                               PB.system.temporal_cache.data[0][2])
PD.write()
pm.IO.collision_list_to_polydata(PB.collisions(), 'collisions%s.vtp'%NAME)


