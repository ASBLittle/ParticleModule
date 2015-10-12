"""Example of particles in double gyre hitting wall."""
import particle_model as pm
import numpy

N = 200

X = numpy.random.random((N, 3))
X[:, 2] = 0
V = numpy.zeros((N, 3))
U = numpy.zeros((N, 3))
GP = numpy.zeros((N, 3))

BOUNDARY = pm.IO.BoundaryData('Gyre_boundary.vtu')
TEMP_CACHE = pm.TemporalCache.TemporalCache('gyre')

PB = pm.Particles.ParticleBucket(X, V, 0.004, 1.0e-4, tc=TEMP_CACHE, U=U, GP=GP,
                                 filename='data.dat', boundary=BOUNDARY,
                                 diameter=1e-4)


PB.write()
for _ in range(1000):
    print PB.time
    print PB.vel[:, 0].ravel().min()
    print PB.fluid_vel[:, 0].ravel().min()
    PB.update()
    PB.write()
PB.outfile.flush()


