"""Example of particles from an options file."""
import particle_model as pm
import numpy

N = 5

X = 0.5+0.25*(numpy.random.random((N, 3))-0.5)

V = numpy.zeros((N,3))

reader=pm.Options.OptionsReader('test.pmml')

def vel(pos):
    """Fluid velocity"""
    return numpy.array((0.1, 0, 0))

def pres(pos):
    """Fluid pressure"""
    return 1.0e3*9.81*(0.1-pos[1])

MESH=pm.IO.get_mesh_from_reader(reader)
pm.IO.make_unstructured_grid(MESH, vel, pres, 0.0, reader.get_name()+'_0.vtu')
pm.IO.make_unstructured_grid(MESH, vel, pres, 5.0, reader.get_name()+'_1.vtu')
SYSTEM=pm.System.get_system_from_reader(reader)
PAR=pm.ParticleBase.get_parameters_from_reader(reader)[0]



PB = pm.Particles.ParticleBucket(X, V, 0.0, 1.0e-2,
                                 system=SYSTEM,
                                 parameters=PAR)


pm.IO.write_level_to_polydata(PB, 0, reader.get_name(), do_average=False)
for i in range(5):
    PB.run(0.1*(i+1),write=False)
    pm.IO.write_level_to_polydata(PB, i+1, reader.get_name(), do_average=False)
