import particle_model as pm
import numpy

N=200

X=numpy.random.random((N,3))
X[:,2]=0
V=numpy.zeros((N,3))
U=numpy.zeros((N,3))
GP=numpy.zeros((N,3))

bnd=pm.IO.boundaryData('Gyre_boundary.vtu')
tc=pm.TemporalCache.TemporalCache('gyre')

pb=pm.Particles.particle_bucket(X,V,0.004,1.0e-4,tc=tc,U=U,GP=GP,
                                filename='data.dat',boundary=bnd,d=1e-4)


pb.write()
for k in range(1000):
    print pb.t
    print pb.v[:,0].ravel().min()
    print pb.u[:,0].ravel().min()
    pb.update()
    pb.write()
pb.file.flush()


