"""Module describing the basic fluids system in which the particles are embedded."""

from particle_model import TemporalCache
from particle_model import Options
from particle_model import IO

from numpy import zeros, empty, bool
from numpy.linalg import norm
import vtk

class System(object):
    """ Class decribes the fixed properties of the underlying system and its
    fluid dynamical solution."""

    def __init__(self, boundary=None, temporal_cache=None, base_name=None, options=None, block=None, velocity_name='Velocity', **kwargs):
        """ Initialise the system class. """
        self.boundary = boundary
        if temporal_cache:
            self.temporal_cache = temporal_cache
        elif base_name:
            self.temporal_cache = TemporalCache.TemporalCache(base_name)
        elif block:
            self.temporal_cache = TemporalCache.FluidityCache(*block,
                                                               velocity_name=velocity_name)
        else:
            self.temporal_cache = None
        self.options = options

        self.gravity = kwargs.get('gravity', zeros(3))
        self.rho = kwargs.get('rho', 1.0e3)
        self.omega = kwargs.get('omega', zeros(3))
        self.viscosity = kwargs.get('viscosity', 1.0e-3)
        self.coeff = kwargs.get('coeff', 0.99)

    def get_rossby(self, velocity, length):
        """ Get the rossby number of the system for a given length scale
        and velocity."""
        return norm(velocity)/(length * 2.0 * norm(self.omega))

    def get_reynolds(self, velocity, length):
        """ Get the Reynolds number of the system for a given lenght scale
        and velocity."""
        return self.rho*norm(velocity)*length/ self.viscosity

    def coefficient_of_restitution(self, particle, cell=None):
        """ Get coefficent of restitution."""

        del particle
        del cell

        return self.coeff

    def in_system(self, points, time):
        """ Check that the points of X are inside the system data """

        out = empty(points.shape[0],bool)

        if self.temporal_cache is None:
            out[:] = True
            return out

        obj = self.temporal_cache(time)[0][0][2]
        loc=vtk.vtkCellLocator()

        if obj.IsA('vtkUnstructuredGrid'):
            loc.SetDataSet(obj)
        else:
            loc.SetDataSet(obj.GetBlock(0))
        loc.BuildLocator()



        for k, point in enumerate(points):
            out[k] = loc.FindCell(point)> -1

        return out
            

def get_system_from_options(options_file=None, boundary_grid=None,
                            block=None, velocity_name='Velocity',dist=None):

    reader=Options.OptionsReader(options_file)

    if boundary_grid is None:

        mesh=IO.GmshMesh()
        mesh.read(reader.get_mesh_filename())

        boundary_grid = IO.make_boundary_from_msh(mesh)

    boundary = IO.BoundaryData(bnd=boundary_grid,
                               outlet_ids=reader.get_outlet_ids(),
                               inlets=reader.get_inlets(), dist=dist)

    if block is None:
        system = System(boundary,base_name=reader.get_name(),
                         gravity=reader.get_gravity(),
                        omega=reader.get_rotation()[0],
                        velocity_name=velocity_name)
    else:
        system = System(boundary,block=block,
                        gravity=reader.get_gravity(),
                        omega=reader.get_rotation()[0],
                        velocity_name=velocity_name)

    return system

def get_system_from_reader(reader, boundary_grid=None):

    if boundary_grid is None:

        mesh=IO.GmshMesh()
        mesh.read(reader.get_mesh_filename())

        boundary_grid = IO.make_boundary_from_msh(mesh)


    boundary = IO.BoundaryData(bnd=boundary_grid)
    system = System(boundary,base_name=reader.get_name(),
                    gravity=reader.get_gravity(),
                    omega=reader.get_rotation()[0])

    return system

    
