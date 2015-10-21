"""Example of particles in double gyre hitting wall."""
import particle_model as pm
import numpy

import pylab as p

S = '160'

N = 20

X = 0.5*(numpy.random.random((N, 3))-0.5)
X[:, 2] += 1.0

V = numpy.zeros((N, 3))
U = numpy.zeros((N, 3))
GP = numpy.zeros((N, 3))

NAME = 'cylinder'

BOUNDARY = pm.IO.BoundaryData('cylinder_boundary.vtu')
TEMP_CACHE = pm.TemporalCache.TemporalCache(NAME)
SYSTEM = pm.System.System(BOUNDARY, gravity=numpy.array((0, 0, -1)))
PAR = pm.Particles.PhysicalParticle(diameter=400e-6)

PB = pm.Particles.ParticleBucket(X, V, 0.0, 5.0e-3, temporal_cache=TEMP_CACHE,
                                 U=U, GP=GP, system=SYSTEM,
                                 e=0.99, parameters=PAR)

TEMP_CACHE.data[1][0] = 100.0

PD = pm.IO.PolyData(NAME+'.vtp')
PD.append_data(PB)
pm.IO.write_level_to_polydata(PB, 0, NAME)

DATA = []
TIME = []

for i in range(30):
    DATA.append(-PB.vel[0, 2])
    TIME.append(PB.time)
    print PB.time
    print 'min, max: pos_z', PB.pos[:, 2].ravel().min(), PB.pos[:, 2].ravel().max()
    print 'min, max: vel_z', PB.vel[:, 2].ravel().min(), PB.vel[:, 2].ravel().max()
    PB.update()
    PD.append_data(PB)
    pm.IO.write_level_to_polydata(PB, i+1, NAME)
DATA.append(-PB.vel[0, 2])
TIME.append(PB.time)

p.plot(TIME, DATA, lw=2)
p.show()

PD.write()
pm.IO.collision_list_to_polydata(PB.collisions(), 'collisions%s.vtp'%NAME)


