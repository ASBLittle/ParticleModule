""" Module containing input-output routines between the particle model and
the file system. Mostly vtk."""

import os
import os.path
import glob
import copy

from particle_model.Debug import profile, logger
from particle_model import Collision
from particle_model import Parallel
from particle_model import vtk_extras

from particle_model.GmshIO import GmshMesh
from particle_model.SolidInteractions import *
from vtk.util import numpy_support
import numpy
from scipy.interpolate import griddata

TYPES_3D = [vtk.VTK_TETRA, vtk.VTK_QUADRATIC_TETRA]
TYPES_2D = [vtk.VTK_TRIANGLE, vtk.VTK_QUADRATIC_TRIANGLE]
TYPES_1D = [vtk.VTK_LINE]

scalar_ugrid_no = {}
vector_ugrid_no = {}
tensor_ugrid_no = {}

ARGV = [0.0, 0.0, 0.0]
WEIGHTS = numpy.zeros(10)
SUB_ID = vtk.mutable(0)

TYPE_DICT = {1 : vtk.VTK_LINE, 2 : vtk.VTK_TRIANGLE, 4 : vtk.VTK_TETRA,
             15 : vtk.VTK_PIXEL}

WRITER = {vtk.VTK_UNSTRUCTURED_GRID:(vtk.vtkXMLPUnstructuredGridWriter
                                     if Parallel.is_parallel()
                                     else vtk.vtkXMLUnstructuredGridWriter),
          vtk.VTK_POLY_DATA:(vtk.vtkXMLPPolyDataWriter
                             if Parallel.is_parallel()
                             else vtk.vtkXMLPolyDataWriter),}

class PolyData(object):
    """ Class storing a living vtkPolyData construction"""

    def __init__(self, filename, fields=None):
        """ Initialize the PolyData instance"""


        if filename.rsplit('.',1)[-1] in ('pvtp','vtp'):
            self.filename = filename
        else:
            if Parallel.is_parallel():
                self.filename = filename +'.pvtp'
            else:
                self.filename = filename +'.vtp'
        self.cell_ids = {}
        self.poly_data = vtk.vtkPolyData()
        self.pnts = vtk.vtkPoints()
        self.pnts.Allocate(0)

        self.fields = fields or {}
        self.fields["Time"] = 1

        self.arrays = {}
        for name, num_comps in self.fields.items():
            array = vtk.vtkDoubleArray()
            array.SetName(name)
            array.SetNumberOfComponents(num_comps)
            self.poly_data.GetPointData().AddArray(array)

    def append_data(self, bucket):
        """ Add the data from the current time level"""

        for particle in bucket:
            ids = self.cell_ids.setdefault(particle, vtk.vtkIdList())

            part_id = self.pnts.InsertNextPoint(particle.pos)

            ids.InsertNextId(part_id)

            for name, num_comps in self.fields.items():
                if name=="Time":
                    self.poly_data.GetPointData().GetArray('Time').InsertNextValue(bucket.time)
                elif name in particle.fields:
                    self.poly_data.GetPointData().GetArray(name).InsertNextValue(particle.fields[name])
                else:
                    data=[]
                    for _ in range(num_comps):
                        data.append(vtk.vtkMath.Nan())
                    self.poly_data.GetPointData().GetArray(name).InsertNextValue(*data)

    def write(self):
        """ Write the staged vtkPolyData to a file."""

        self.poly_data.SetPoints(self.pnts)

        self.poly_data.Allocate(len(self.cell_ids))
        for cell_id in self.cell_ids.values():
            self.poly_data.InsertNextCell(vtk.VTK_LINE, cell_id)

        writer = WRITER[vtk.VTK_POLY_DATA]()
        writer.SetFileName(self.filename)
        if Parallel.is_parallel():
            writer.SetNumberOfPieces(Parallel.get_size())
            writer.SetStartPiece(Parallel.get_rank())
            writer.SetEndPiece(Parallel.get_rank())
            if vtk.vtkVersion.GetVTKMajorVersion()==6:
                writer.SetWriteSummaryFile(Parallel.get_rank()==0)
            else:
                controller = vtk.vtkMPIController()
                controller.SetCommunicator(vtk.vtkMPICommunicator.GetWorldCommunicator())
                writer.SetController(controller)
        if vtk.vtkVersion.GetVTKMajorVersion() < 6:
            writer.SetInput(self.poly_data)
        else:
            writer.SetInputData(self.poly_data)
        writer.Write()

        if Parallel.is_parallel():
            make_subdirectory(self.filename)

class BoundaryData(object):
    """ Class storing the boundary data for the problem"""
    def __init__(self, filename=None, bnd=None, outlet_ids=[], inlets=[], dist=None):
        """Class containing the information about the boundary of the domain.

        Args:
            filename (str): Name of the file containing the
            vtkUnstructuredGrid denoting the boundary of the domain."""

        self.reader = vtk.vtkXMLUnstructuredGridReader()
        self.bndl = vtk.vtkCellLocator()
        self.geom_filter = vtk.vtkGeometryFilter()
        self.outlet_ids = outlet_ids
        self.inlets = inlets
        self.dist = dist

        open_ids = [] + self.outlet_ids
        for inlet in self.inlets:
            open_ids += inlet.surface_ids

        if filename is not None:
            self.update_boundary_file(filename, open_ids)
        else:
            self.bnd = bnd
            if self.dist:
                self.phys_bnd = move_boundary_through_normal(self.bnd, self.dist,
                                                             Ids=open_ids)
            else:
                self.phys_bnd = self.bnd

            if vtk.vtkVersion.GetVTKMajorVersion() < 6:
                self.geom_filter.SetInput(self.phys_bnd)
            else:
                self.geom_filter.SetInputData(self.phys_bnd)
            self.geom_filter.Update()

            self.bndl.SetDataSet(self.geom_filter.GetOutput())
            self.bndl.BuildLocator()

    def update(self, boundary):
        """ Update the boundary data from a new object."""

        if boundary.IsA('vtkMultiBlockDataSet'):
            self.bnd = boundary.GetBlock(boundary.GetNumberOfBlocks()-1)
        else:
            self.bnd = boundary

        self.bndl.SetDataSet(self.bnd)
        self.bndl.BuildLocator()

        self.bndl.SetDataSet(self.bnd)
        self.bndl.BuildLocator()

    def update_boundary_file(self, infile, open_ids=None):
        """ Update the boundary data from the file."""

        open_ids = open_ids or []

        if isinstance(infile, str):
            if not os.path.isfile(infile):
                logger.error(os.getcwd())
                raise OSError
            self.reader.SetFileName(infile)
            self.reader.Update()
            self.bnd = self.reader.GetOutput()
        else:
            self.bnd = infile

        if self.dist:
            self.phys_bnd = move_boundary_through_normal(self.bnd, self.dist,
                                                         Ids=open_ids)
        else:
            self.phys_bnd = self.bnd

        if vtk.vtkVersion.GetVTKMajorVersion() < 6:
            self.geom_filter.SetInput(self.phys_bnd)
        else:
            self.geom_filter.SetInputData(self.phys_bnd)
        self.geom_filter.Update()

        self.bndl.SetDataSet(self.geom_filter.GetOutput())
        self.bndl.BuildLocator()

    def rebuild_locator(self):
        """ Rebuild the locator information"""
        self.bndl.BuildLocator()



def clean_unstructured_grid(ugrid):
    """Collapse a vtu produced from a discontinuous grid back down to the continuous space.

    Args:
    ugrid (vtkUnstructuredGrid): the input discontinuous grid

    Results
    out_grid (vtkUnstructuredGrid): A continuous grid"""

    merge_points = vtk.vtkMergePoints()
    out_grid = vtk.vtkUnstructuredGrid()

    for i in range(ugrid.GetNumberOfPoints()):
        merge_points.InsertUniquePoint(ugrid.GetPoints().GetPoint(i))

    merge_points.BuildLocator()

    pts = vtk.vtkPoints()
    pts.DeepCopy(merge_points.GetPoints())
    out_grid.SetPoints(pts)

    for i in range(ugrid.GetNumberOfCells()):
        cell = ugrid.GetCell(i)
        cell_ids = cell.GetPointIds()

        for j in range(cell.GetNumberOfPoints()):

            original_point = cell.GetPoints().GetPoint(j)
            cell_ids.SetId(j,
                           merge_points.FindClosestInsertedPoint(original_point))

        out_grid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())


    out_grid.GetCellData().DeepCopy(ugrid.GetCellData())

    return out_grid

def extract_boundary(ugrid):
    """Extract the boundary elements from an unstructured grid, provided it already contains them.

    Args:

    ugrid (vtkUnstructuredGrid): The grid with which to work.

    Results:

    out_grid (vtkUnstructuredGrid): Grid containing the boundary of ugrid"""

    out_grid = vtk.vtkUnstructuredGrid()
    pts = vtk.vtkPoints()
    pts.DeepCopy(ugrid.GetPoints())
    out_grid.SetPoints(pts)
    out_grid.GetCellData().CopyStructure(ugrid.GetCellData())

    celltypes = vtk.vtkCellTypes()

    ugrid.GetCellTypes(celltypes)

    if any([celltypes.IsType(ct) for ct in TYPES_3D]):
        dim = 3
    elif any([celltypes.IsType(ct) for ct in TYPES_2D]):
        dim = 2
    elif any([celltypes.IsType(ct) for ct in TYPES_1D]):
        dim = 1
    else:
        dim = 0

    def get_dimension(cell):
        """Get dimensionality of cell."""
        cell_type = cell.GetCellType()
        if cell_type in TYPES_3D:
            return 3
        if cell_type in TYPES_2D:
            return 2
        if cell_type in TYPES_1D:
            return 1
        #otherwise
        return 0

    ncells = ugrid.GetNumberOfCells()
    ncda = ugrid.GetCellData().GetNumberOfArrays()

    for i in range(ncda):
        out_grid.GetCellData().GetArray(i).SetName(ugrid.GetCellData().GetArray(i).GetName())

    cell_data = ugrid.GetCellData()
    for i in range(ncells):
        cell = ugrid.GetCell(i)
        if dim > get_dimension(cell):
            out_grid.InsertNextCell(cell.GetCellType(),
                                    cell.GetPointIds())
            for j in range(ncda):
                out_data = out_grid.GetCellData().GetArray(j)
                out_data.InsertNextTuple(cell_data.GetArray(j).GetTuple(i))

    return out_grid

def write_bucket_to_polydata(bucket):

    """ Output the points of a bucket to a vtkPoints object. """

    poly_data = vtk.vtkPolyData()
    pnts = vtk.vtkPoints()
    pnts.Allocate(0)
    poly_data.SetPoints(pnts)
    poly_data.Allocate(len(bucket))

    for positions, in bucket.pos():
        pixel = vtk.vtkPixel()
        pixel.GetPointIds().InsertId(0,
                                     poly_data.GetPoints().InsertNextPoint(*positions))
        poly_data.InsertNextCell(pixel.GetCellType(), pixel.GetPointIds())

def write_bucket_to_points(bucket):

    """ Output the basic data of a bucket to a vtkPolyData object. """

    pnts = vtk.vtkPoints()
    pnts.Allocate(len(bucket))

    for positions in bucket.pos():
        pnts.InsertNextPoint(*positions)

    return pnts

def write_level_to_polydata(bucket, level, basename=None, do_average=False,
                            field_data=None, **kwargs):

    """Output a time level of a particle bucket to a vtkPolyData (.vtp) files.

    Each file contains one time level of the data, and are numbered sequentially.
    Within each file, each particle is written to seperate pixel.

    Args:
         bucket   (ParticleBucket):
        level    (int):
        basename (str): String used in the construction of the file series.
        The formula is of the form basename_0.vtp, basename_1.vtp,..."""

    del kwargs
    field_data = field_data or {}

    poly_data = vtk.vtkPolyData()
    pnts = vtk.vtkPoints()
    pnts.Allocate(0)
    poly_data.SetPoints(pnts)
    poly_data.Allocate(len(bucket))

    outtime = vtk.vtkDoubleArray()
    outtime.SetName('Time')
    outtime.Allocate(1)

    particle_id = vtk.vtkDoubleArray()
    particle_id.SetName('ParticleID')
    particle_id.Allocate(len(bucket))

    live = vtk.vtkDoubleArray()
    live.SetName('Live')
    live.Allocate(len(bucket))

    plive = bucket.system.in_system(bucket.pos(), len(bucket), bucket.time)

    for _, par in zip(plive, bucket):
        particle_id.InsertNextValue(hash(par))
        if _:
            live.InsertNextValue(1.0)
        else:
            live.InsertNextValue(0.0)


    velocity = vtk.vtkDoubleArray()
    velocity.SetNumberOfComponents(3)
    velocity.Allocate(len(bucket))
    velocity.SetName('Particle Velocity')

    for positions, vel in zip(bucket.pos(), bucket.vel()):
        pixel = vtk.vtkPixel()
        pixel.GetPointIds().InsertId(0,
                                     poly_data.GetPoints().InsertNextPoint(positions[0],
                                                                           positions[1],
                                                                           positions[2]))
        velocity.InsertNextTuple3(vel[0], vel[1], vel[2])
        poly_data.InsertNextCell(pixel.GetCellType(), pixel.GetPointIds())

    outtime.InsertNextValue(bucket.time)

    poly_data.GetFieldData().AddArray(outtime)
    poly_data.GetPointData().AddArray(velocity)
    poly_data.GetPointData().AddArray(particle_id)
    poly_data.GetCellData().AddArray(live)

    for name, num_comps in field_data.items():
        _ = vtk.vtkDoubleArray()
        _.SetName(name)
        _.SetNumberOfComponents(num_comps)
        _.Allocate(len(bucket))
        for k, particle in enumerate(bucket):
            _.InsertNextValue(particle.fields[name])
        poly_data.GetPointData().AddArray(_)

    if do_average:
        gsp = calculate_averaged_properties_cpp(poly_data)

    if Parallel.is_parallel():
        file_ext = 'pvtp'
    else:
        file_ext = 'vtp'

    write_to_file(poly_data, "%s_%d.%s"%(basename, level, file_ext))

    if do_average:
        return gsp

def write_level_to_ugrid(bucket, level, basename, model, **kwargs):
    """ Output a time level of a bucket to a vtkXMLUnstructuredGrid (.vtu) file.

    Particle volumes are averaged over the cells in the model using a control volume
    approach."""

    del kwargs

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.DeepCopy(model)

#    mass = numpy.zeros(ugrid.GetNumberOfPoints())
#    for _ in range(ugrid.GetNumberOfCells()):
#        cell = ugrid.GetCell(_)
#        dummy = cell.GetPoints().GetPoint
#        dummy_id = cell.GetPointIds().GetId
#
#        measure = abs(cell.TriangleArea(dummy(0), dummy(1), dummy(2)))
#
#        for k in range(cell.GetNumberOfPoints()):
#            mass[dummy_id(k)] += measure/cell.GetNumberOfPoints()

    locator = vtk.vtkPointLocator()
    locator.SetDataSet(ugrid)
    locator.BuildLocator()

    volume = numpy.zeros(ugrid.GetNumberOfPoints())
    temperature = numpy.zeros(ugrid.GetNumberOfPoints())
    solid_pressure = numpy.zeros(ugrid.GetNumberOfPoints())
    velocity = numpy.zeros((ugrid.GetNumberOfPoints(), 3))

    LENGTH = 0.1
    MULTIPLIER = 1.e3

    for dummy_particle in bucket:
        point_list = vtk.vtkIdList()
        locator.FindPointsWithinRadius(LENGTH, dummy_particle.pos, point_list)

        for _ in range(point_list.GetNumberOfIds()):
            point_index = point_list.GetId(_)

            rad2 = numpy.sum((numpy.array(ugrid.GetPoints().GetPoint(point_index))
                              -dummy_particle.pos)**2)
            rad2 /= LENGTH

            volume[point_index] += (1.0/6.0*numpy.pi*dummy_particle.parameters.diameter**3
                                    *numpy.exp(-rad2**2)*MULTIPLIER)
            velocity[point_index, :] += (dummy_particle.vel*1.0/6.0*numpy.pi
                                         *dummy_particle.parameters.diameter**3
                                         *numpy.exp(-rad2**2)*MULTIPLIER)

    volume /= 0.5*LENGTH**2*(1.0-numpy.exp(-1.0**2))
    velocity /= 0.5*LENGTH**2*(1.0-numpy.exp(-1.0**2))

    for dummy_particle in bucket:
        point_list = vtk.vtkIdList()
        locator.FindPointsWithinRadius(LENGTH, dummy_particle.pos, point_list)

        rad2 = numpy.sum((numpy.array(ugrid.GetPoints().GetPoint(point_index))
                          -dummy_particle.pos)**2)
        rad2 /= LENGTH

        for _ in range(point_list.GetNumberOfIds()):
            point_index = point_list.GetId(_)
            c = numpy.sum((dummy_particle.vel-velocity[point_index, :])**2)

            temperature[point_index] += (c*1.0/6.0*numpy.pi
                                         *dummy_particle.parameters.diameter**3
                                         *numpy.exp(-rad2**2)*MULTIPLIER)

    for _ in range(ugrid.GetNumberOfPoints()):

        solid_pressure[_] = (dummy_particle.parameters.rho*volume[_]
                             *radial_distribution_function(volume[_])*temperature[_])

    data = [vtk.vtkDoubleArray()]
    data[0].SetName('SolidVolumeFraction')
    data.append(vtk.vtkDoubleArray())
    data[1].SetName('SolidVolumeVelocity')
    data[1].SetNumberOfComponents(3)
    data.append(vtk.vtkDoubleArray())
    data[2].SetName('GranularTemperature')
    data.append(vtk.vtkDoubleArray())
    data[3].SetName('SolidPressure')

    for _ in range(ugrid.GetNumberOfPoints()):
        data[0].InsertNextValue(volume[_])
        if volume[_] > 1.0e-20:
            data[1].InsertNextTuple3(*(velocity[_]/volume[_]))
        else:
            data[1].InsertNextTuple3(*(0*velocity[_]))
        data[2].InsertNextValue(temperature[_])
        data[3].InsertNextValue(solid_pressure[_])
        ugrid.GetPointData().GetScalars('Time').SetValue(_, bucket.time)

    for _ in data:
        ugrid.GetPointData().AddArray(_)

    write_to_file(ugrid, "%s_out_%d.vtu"%(basename, level))

def update_collision_polydata(bucket, base_name, **kwargs):
    """ Update collisions from data in the particle bucket."""

    if Parallel.is_parallel():
        fext = 'pvtp'
    else:
        fext = 'vtp'

    collision_list_to_polydata(bucket.collisions(), base_name+'_collisions.'+fext)

def collision_list_to_polydata(col_list, outfile,
                               model=Collision.mclaury_mass_coeff, **kwargs):
    """Convert collision data to a single vtkPolyData (.vtp) files.

    Each particle is written to seperate cell.

    Args:
        outfile (str):  Filename of the output PolyDataFile. The extension .vtp
        is NOT added automatically."""

    poly_data = vtk.vtkPolyData()
    pnts = vtk.vtkPoints()
    pnts.Allocate(0)
    poly_data.SetPoints(pnts)
    poly_data.Allocate(len(col_list))

    time = vtk.vtkDoubleArray()
    time.SetName('Time')
    wear = vtk.vtkDoubleArray()
    wear.SetName('Wear')
    normal = vtk.vtkDoubleArray()
    normal.SetNumberOfComponents(3)
    normal.SetName('Normal')

    for col in col_list:
        pixel = vtk.vtkPixel()
        pixel.GetPointIds().InsertId(0,
                                     poly_data.GetPoints().InsertNextPoint(col.pos[0],
                                                                           col.pos[1],
                                                                           col.pos[2]))
        time.InsertNextValue(col.time)
        wear.InsertNextValue(model(col, **kwargs))
        normal.InsertNextTuple3(*col.normal)
        poly_data.InsertNextCell(pixel.GetCellType(), pixel.GetPointIds())

    poly_data.GetPointData().AddArray(time)
    poly_data.GetPointData().AddArray(wear)
    poly_data.GetPointData().AddArray(normal)

    write_to_file(poly_data, outfile)


def get_linear_cell(cell):
    """ Get equivalent linear cell to vtkCell cell"""
    if cell.GetCellType() in (vtk.VTK_POLY_LINE,):
        linear_cell = vtk.vtkLine()
        linear_cell.GetPoints().SetPoint(0, cell.GetPoints().GetPoint(0))
        linear_cell.GetPoints().SetPoint(1, cell.GetPoints().GetPoint(1))
    elif cell.GetCellType() in (vtk.VTK_QUADRATIC_TRIANGLE,):
        linear_cell = vtk.vtkTriangle()
        linear_cell.GetPoints().SetPoint(0, cell.GetPoints().GetPoint(0))
        linear_cell.GetPoints().SetPoint(1, cell.GetPoints().GetPoint(1))
        linear_cell.GetPoints().SetPoint(2, cell.GetPoints().GetPoint(2))
    elif cell.GetCellType() in (vtk.VTK_QUADRATIC_TETRA,):
        linear_cell = vtk.vtkTetra()
        linear_cell.GetPoints().SetPoint(0, cell.GetPoints().GetPoint(0))
        linear_cell.GetPoints().SetPoint(1, cell.GetPoints().GetPoint(1))
        linear_cell.GetPoints().SetPoint(2, cell.GetPoints().GetPoint(2))
        linear_cell.GetPoints().SetPoint(3, cell.GetPoints().GetPoint(3))
    else:
        linear_cell = cell

    return linear_cell

def test_in_cell(cell, position):
    """ Check if point is in vtk cell"""

    linear_cell = get_linear_cell(cell)
    dim = linear_cell.GetNumberOfPoints()-1
    ppos = numpy.zeros(linear_cell.GetNumberOfPoints())
    dummy_func = cell.GetPoints().GetPoint
    args = [dummy_func(i)[:dim] for i in range(1, dim+1)]
    args.append(dummy_func(0)[:dim])
    args.append(ppos)
    linear_cell.BarycentricCoords(position[:dim], *args)

    out = cell.GetParametricDistance(ppos[:3])
    if out > 0:
        logger.warn('point %s'%position)
        logger.warn('outside cell by %s'%out)
    return out == 0

def write_to_file(vtk_data, outfile):
    """ Wrapper around the various VTK writer routines"""

    writer = WRITER[vtk_data.GetDataObjectType()]()
    writer.SetFileName(outfile)
    if Parallel.is_parallel():
        writer.SetNumberOfPieces(Parallel.get_size())
        writer.SetStartPiece(Parallel.get_rank())
        writer.SetEndPiece(Parallel.get_rank())
        if vtk.vtkVersion.GetVTKMajorVersion()==6:
            writer.SetWriteSummaryFile(Parallel.get_rank()==0)
        else:
            controller = vtk.vtkMPIController()
            controller.SetCommunicator(vtk.vtkMPICommunicator.GetWorldCommunicator())
            writer.SetController(controller)
    if vtk.vtkVersion.GetVTKMajorVersion()<6:
        writer.SetInput(vtk_data)
    else:
        writer.SetInputData(vtk_data)        
    writer.Write()

    if Parallel.is_parallel():
        make_subdirectory(outfile)

def make_subdirectory(fname):
    Parallel.barrier()
    if Parallel.get_rank()==0:
        base_name, file_type = fname.rsplit('.',1)
        if not os.path.isdir(base_name):
            os.mkdir(base_name)
        for _ in glob.glob(base_name+'_*.%s'%file_type[1:]):
            os.rename(_, base_name+'/'+_)
            with open(fname) as summary_file:
                new_text = summary_file.read().replace(_, base_name+'/'+_)
            with open(fname,'w') as summary_file:
                summary_file.write(new_text)

def make_unstructured_grid(mesh, velocity, pressure, time, outfile=None):
    """Given a mesh (in Gmsh format), velocity and pressure fields, and a
    time level, store the data in a vtkUnstructuredGridFormat."""

    pnts = vtk.vtkPoints()
    pnts.Allocate(len(mesh.nodes))

    node2id = {}
    for k, point in mesh.nodes.items():
        node2id[k] = pnts.InsertNextPoint(point)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(pnts)

    for element in mesh.elements.values():
        id_list = vtk.vtkIdList()
        for node in element[2]:
            id_list.InsertNextId(node2id[node])

        ugrid.InsertNextCell(TYPE_DICT[element[0]], id_list)

    data = []
    data.append(vtk.vtkDoubleArray())
    data[0].SetNumberOfComponents(3)
    data[0].Allocate(3*pnts.GetNumberOfPoints())
    data[0].SetName('Velocity')

    data.append(vtk.vtkDoubleArray())
    data[1].Allocate(pnts.GetNumberOfPoints())
    data[1].SetName('Pressure')

    data.append(vtk.vtkDoubleArray())
    data[2].Allocate(pnts.GetNumberOfPoints())
    data[2].SetName('Time')

    for k in range(len(mesh.nodes)):
        if hasattr(velocity, '__call__'):
            data[0].InsertNextTuple3(*velocity(ugrid.GetPoints().GetPoint(k)))
        else:
            data[0].InsertNextTuple3(*velocity[k, :])
        if hasattr(pressure, '__call__'):
            data[1].InsertNextValue(pressure(ugrid.GetPoints().GetPoint(k)))
        else:
            data[1].InsertNextValue(pressure[k])
        data[2].InsertNextValue(time)


    for _ in data:
        ugrid.GetPointData().AddArray(_)

    if outfile:
        write_to_file(ugrid, outfile)

    return ugrid

def get_boundary_from_block(mblock):
    """ Get properly formed boundary data from a fluidity vtk block."""
    
    ugrid = None

    for _ in range(mblock.GetNumberOfBlocks()):
        if mblock.GetMetaData(_).Get(vtk.vtkCompositeDataSet.NAME()) == 'Boundary':
            ugrid = mblock.GetBlock(_)

    if not ugrid:
        return None

    vgrid = vtk.vtkUnstructuredGrid()

    pts = vtk.vtkPoints()
    pts.DeepCopy(ugrid.GetPoints())

    vgrid.SetPoints(pts)

    surface_ids=vtk.vtkIntArray()
    surface_ids.SetName('SurfaceIds')
    surface_ids.SetNumberOfComponents(1)
    surface_ids.SetNumberOfTuples(ugrid.GetNumberOfCells())
    
    sids = ugrid.GetCellData().GetArray('SurfaceIds')

    for _ in range(ugrid.GetNumberOfCells()):
        cell = ugrid.GetCell(_)
        vgrid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())
        surface_ids.SetValue(_, int(sids.GetValue(_)))

    vgrid.GetCellData().AddArray(surface_ids)

    return vgrid

def make_boundary_from_msh(mesh, outfile=None):

    """Given a mesh (in Gmsh format), store the boundary data in a vtkUnstructuredGridFormat."""

    pnts = vtk.vtkPoints()
    pnts.Allocate(len(mesh.nodes))

    node2id = {}
    for k, point in mesh.nodes.items():
        node2id[k] = pnts.InsertNextPoint(point)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(pnts)

    surface_ids=vtk.vtkIntArray()
    surface_ids.SetName('SurfaceIds')

    for element in mesh.elements.values():
        id_list = vtk.vtkIdList()
        for node in element[2]:
            id_list.InsertNextId(node2id[node])
        ugrid.InsertNextCell(TYPE_DICT[element[0]], id_list)
        if element[1]:
            surface_ids.InsertNextValue(element[1][0])


    ugrid.GetCellData().AddArray(surface_ids)

    ugrid_bnd = extract_boundary(ugrid)

    if outfile:
        write_to_file(ugrid_bnd, outfile)

    return ugrid_bnd

def interpolate_collision_data(col_list, ugrid, method='nearest'):
    """ Interpolate the wear data from col_list onto the surface described
    in the vtkUnstructuredGrid ugrid """

    out_pts = numpy.array([ugrid.GetPoint(i) for i in range(ugrid.GetNumberOfPoints())])

    data_pts = numpy.array([c.pos for c in col_list])
    wear_pts = numpy.array([c.get_wear() for c in col_list])

    wear_on_grid = griddata(data_pts, wear_pts, out_pts, method=method)
#    rbfi = Rbf(data_pts[:,0], data_pts[:, 1], data_pts[:, 2], wear_pts,
#    function='gaussian', epsilon=0.1,smooth=3.0)
#    wear_on_grid = rbfi(out_pts[:,0], out_pts[:, 1], out_pts[:, 2])

    wear_vtk = vtk.vtkDoubleArray()
    wear_vtk.SetName('Wear')

    for wear in wear_on_grid:
        wear_vtk.InsertNextValue(wear)

    ugrid.GetPointData().AddArray(wear_vtk)


def get_linear_block(infile):
    if infile.IsA('vtkUnstructuredGrid'):
        return infile
    elif infile.IsA('vtkMultiBlockDataSet'):
        return infile.GetBlock(0)
    else:
        raise AttributeError

@profile
def get_tensor(infile, name, index, pcoords):

    global tensor_ugrid_no

    if infile.IsA('vtkUnstructuredGrid'):
        ids = infile.GetCell(index).GetPointIds()
        data = infile.GetPointData().GetVectors(name)
        ugrid = infile
    else:
        if tensor_ugrid_no.setdefault(name, None):
            ids = infile.GetBlock(tensor_ugrid_no[name]).GetCell(index).GetPointIds()
            data = infile.GetBlock(tensor_ugrid_no[name]).GetPointData().GetScalars(name)
            ugrid = infile.GetBlock(tensor_ugrid_no[name])
        else:
            for _ in range(infile.GetNumberOfBlocks()):
                if infile.GetBlock(_).GetPointData().HasArray(name):
                    ids = infile.GetBlock(_).GetCell(index).GetPointIds()
                    data = infile.GetBlock(_).GetPointData().GetVectors(name)
                    ugrid = infile.GetBlock(_)
                    tensor_ugrid_no[name] = _
                    break

#    ugrid.GetCell(index).EvaluateLocation(SUB_ID, pcoords, ARGV, WEIGHTS)
    ugrid.GetCell(index).InterpolateFunctions(pcoords, WEIGHTS[:ids.GetNumberOfIds()])

    dim = numpy.sqrt(data.GetNumberOfComponents())
    out = numpy.zeros((dim,dim), float)
    for _ in range(ids.GetNumberOfIds()):
        out[:] += WEIGHTS[_]*numpy.array(data.GetTuple(ids.GetId(_))).reshape((dim,dim))
    return out 

@profile
def get_vector(infile, data, name, index, pcoords):

    global vector_ugrid_no

    if infile.IsA('vtkUnstructuredGrid'):
        ids = infile.GetCell(index).GetPointIds()
        if not data:
            data = infile.GetPointData().GetVectors(name)
        ugrid = infile
    else:
        if vector_ugrid_no.setdefault(name,None):
            ids = infile.GetBlock(vector_ugrid_no[name]).GetCell(index).GetPointIds()
            if not data:
                data = infile.GetBlock(vector_ugrid_no[name]).GetPointData().GetScalars(name)
            ugrid = infile.GetBlock(vector_ugrid_no[name])
        else:
            for _ in range(infile.GetNumberOfBlocks()):
                if infile.GetBlock(_).GetPointData().HasArray(name):
                    ids = infile.GetBlock(_).GetCell(index).GetPointIds()
                    if not data:
                        data = infile.GetBlock(_).GetPointData().GetVectors(name)
                    ugrid = infile.GetBlock(_)
                    vector_ugrid_no[name] = _
                    break
            if vector_ugrid_no[name] is None: return None
    if not data: return numpy.zeros(3) 

#    ugrid.GetCell(index).EvaluateLocation(SUB_ID, pcoords, ARGV, WEIGHTS)
    ugrid.GetCell(index).InterpolateFunctions(pcoords, WEIGHTS[:ids.GetNumberOfIds()])

    out = vtk_extras.vInterpolate(data, ids, WEIGHTS) 
    return out 
            
@profile
def get_scalar(infile, data, name, index):
        
    global scalar_ugrid_no

    if infile.IsA('vtkUnstructuredGrid'):
        ids = infile.GetCell(index).GetPointIds()
        if not data:
            data = infile.GetPointData().GetScalars(name)
        scalar_ugrid = infile
    else:
        if scalar_ugrid_no.setdefault(name,None):
            ids = infile.GetBlock(scalar_ugrid_no[name]).GetCell(index).GetPointIds()
            if not data:
                data = infile.GetBlock(scalar_ugrid_no[name]).GetPointData().GetScalars(name)
        else:
            for _ in range(infile.GetNumberOfBlocks()):
                if infile.GetBlock(_).GetPointData().HasArray(name):
                    ids = infile.GetBlock(_).GetCell(index).GetPointIds()
                    if not data:
                        data = infile.GetBlock(_).GetPointData().GetScalars(name)
                    scalar_ugrid_no[name] = _
                    break
    if data:
        out = numpy.empty(ids.GetNumberOfIds(),float)
        for _ in range(ids.GetNumberOfIds()):
            out[_]=data.GetValue(ids.GetId(_))
    else:
        out = None
    return out 

def get_mesh_from_reader(reader):
    mesh = GmshMesh()
    mesh.read(reader.get_mesh_filename())
    return mesh

def get_boundary_from_fluidity_mesh(positions):
    """Temporary method for testing"""

    faces={}

    surface_nodes=set()
    
    cells=[]
    
    for i in range(positions.element_count):
        nodes=positions.ele_nodes(i)
        for node in nodes:
            tmp=copy.copy(nodes)
            tmp.remove(node)
            tmp=tuple(sorted(tmp))
            faces[tmp]=not faces.get(tmp, False)

    for key, val in faces.items():
        if val:
            for node in key:
                surface_nodes.add(node)
            cells.append(key)
        
    ugrid = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    points.Allocate(len(surface_nodes))
    
    nodemap={}
    for key, node in enumerate(surface_nodes):
        X = numpy.empty(3,float)
        X[:positions.dimension] = positions.node_val(node)
        points.InsertNextPoint(X)
        nodemap[node]=key
    ugrid.SetPoints(points)

    for cell in cells:
        if len(cell) == 3:
            cell_type=vtk.VTK_TRIANGLE
        else:
            cell_type=vtk.VTK_LINE

        pntids=vtk.vtkIdList()
        for node in cell:
            pntids.InsertNextId(nodemap[node])
        ugrid.InsertNextCell(cell_type,pntids)

    return ugrid

def move_boundary_through_normal(ugrid, distance, Ids=[]):

    r = numpy.zeros((ugrid.GetNumberOfPoints(), 3))
    area = numpy.zeros(ugrid.GetNumberOfPoints())
    n = numpy.zeros((ugrid.GetNumberOfCells(), 3))
    local_area =  numpy.zeros(ugrid.GetNumberOfCells())


    
    if ugrid.GetCell(0).GetClassName() == 'vtkLine':
    

        for i in range(ugrid.GetNumberOfCells()):

            cell = ugrid.GetCell(i)

            x = numpy.array(cell.GetPoints().GetPoint(1))
            p = x - numpy.array(cell.GetPoints().GetPoint(0))

            local_area = numpy.sqrt(p.dot(p))
            p = p/local_area

            pid0 = cell.GetPointId(0)
            pid1 = cell.GetPointId(1)

            n[i,0] = p[1]
            n[i,1] = -p[0]
            
            if abs(n[i,0])>0.0:
                if x[0]>0.11:
                    n[i,0]=1.0
                else:
                    n[i,0]=-1.0
            else:
                if x[1]>0.0:
                    n[i,1]=1.0
                else:
                    n[i,1]=-1.0
                

    else:

        gf = vtk.vtkGeometryFilter()
        if vtk.vtkVersion.GetVTKMajorVersion()<6:
            gf.SetInput(ugrid)
        else:
            gf.SetInput(ugrid)
        gf.Update()
        
        pd = gf.GetOutput()

        pn = vtk.vtkPolyDataNormals()
        if vtk.vtkVersion.GetVTKMajorVersion()<6:
            pn.SetInput(pd)
        else:
            pn.SetInputData(pd)
        pn.ComputeCellNormalsOn()
        pn.AutoOrientNormalsOn()
        pn.Update()
        norms = pn.GetOutput().GetCellData().GetNormals()
        
        N=10

        for k in range(N):
           
            r0 = numpy.zeros((ugrid.GetNumberOfPoints(), 3))

            physical_ids = ugrid.GetCellData().GetScalars("PhysicalIds")

            for i in range(ugrid.GetNumberOfCells()):

                if physical_ids.GetValue(i) in Ids:
                    continue

                cell = ugrid.GetCell(i)
                n[i,:] = norms.GetTuple(i)

                p0 = (numpy.array(cell.GetPoints().GetPoint(0))-
                      k/float(N-1)*distance*r[cell.GetPointId(0),:])
                p1 = (numpy.array(cell.GetPoints().GetPoint(1))-
                      k/float(N-1)*distance*r[cell.GetPointId(1),:])
                p2 = (numpy.array(cell.GetPoints().GetPoint(2))-
                      k/float(N-1)*distance*r[cell.GetPointId(2),:])
                    
                local_area[i] = cell.TriangleArea(p0,p1,p2)

                for j in range(cell.GetNumberOfPoints()):
               
                    pid = cell.GetPointId(j) 
                    r0[pid,:] += local_area[i] * n[i,:]
                    area[pid] += local_area[i]

            r[:,:] = 0.0

            for i in range(ugrid.GetNumberOfCells()):

                if physical_ids.GetValue(i) in Ids:
                    continue

                cell = ugrid.GetCell(i)

                for j in range(cell.GetNumberOfPoints()):

                    pid = cell.GetPointId(j) 
                    
                    r[pid,:] += local_area[i] * n[i,:] / n[i,:].dot(r0[pid,:])

        for i in range(ugrid.GetNumberOfPoints()):

            x=numpy.empty(3)
            ugrid.GetPoints().GetPoint(i,x)
            ugrid.GetPoints().SetPoint(i,x-distance*r[i,:])

    return ugrid

    
def get_real_x(cell, locx):
    ### Return physical coordinate of cell corresponding to locx in vcell
    ### Note the ordering used here is because vtk is weird.

    if cell.GetCellType() == vtk.VTK_LINE:
        return numpy.array(cell.GetPoint(0))*(1.0-locx[0])+numpy.array(cell.GetPoint(1))*locx[0]
    else:
        return numpy.array(cell.GetPoint(0))*(1.0-locx[0]-locx[1])+numpy.array(cell.GetPoint(1))*locx[0]+numpy.array(cell.GetPoint(2))*locx[1]
