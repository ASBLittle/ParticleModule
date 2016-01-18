""" Module to convert from the state_type to vtk. """

import vtk
from vtk.util import numpy_support
import numpy

try:
    from mpi4py import MPI
except ImportError:
    pass
    

import os
import errno


NUM_DICT = {('lagrangian', 2, 3) : [0,1,2],
             ('lagrangian', 2, 6) : [0,3,1,5,4,2],
             ('lagrangian', 3, 4) : [0,1,2,3],
             ('lagrangian', 3, 10) : [0,4,1,6,5,2,7,8,9,3]}

CELL_DICT ={ ('lagrangian', 1, 2) : vtk.VTK_LINE,
             ('lagrangian', 1, 3) : vtk.VTK_QUADRATIC_EDGE,
             ('lagrangian', 2, 3) : vtk.VTK_TRIANGLE,
             ('lagrangian', 2, 6) : vtk.VTK_QUADRATIC_TRIANGLE,
             ('lagrangian', 3, 4) : vtk.VTK_TETRA,
             ('lagrangian', 3, 10) : vtk.VTK_QUADRATIC_TETRA}

def fluidity_to_mblock(state):

    """Convert a fluidity python state into a vtk multiblock dataset. 

    Individual blocks are created for:

    1. P0DG and P1 continuous fields
    2. P1DG fields (if present)
    3. P2 fields (if present)
    4. P2DG fields (if present)

    Returns the vtkMultiBlockDataSet object."""

    mblock = vtk.vtkMultiBlockDataSet()

    # Deal with the P1 data (always present

    ugrid = fluidity_to_ugrid_p1(state, is_p1, is_p0,
                                 state.vector_fields['Coordinate'])
    mblock.SetBlock(0, ugrid)
    mblock.GetMetaData(0).Set(vtk.vtkCompositeDataSet.NAME(), 'P1CG') 

    dummy = 1

    for name, mesh_test in (('P1DG', is_p1dg),
                            ('P2CG', is_p2),
                            ('P2DG', is_p2dg)):
        ugrid = fluidity_to_ugrid_by_mesh(state, mesh_test)
        if ugrid:
            mblock.SetBlock(dummy, ugrid)
            mblock.GetMetaData( dummy ).Set( vtk.vtkCompositeDataSet.NAME(),
                                             name )
            dummy += 1

    ugrid = fluidity_to_ugrid_p1(state,is_surface_p1, is_surface_p0,
                                 state.vector_fields['SurfaceCoordinate'])
    mblock.SetBlock(dummy, ugrid)
    mblock.GetMetaData(dummy).Set(vtk.vtkCompositeDataSet.NAME(), 'Boundary' )

    writer=vtk.vtkXMLMultiBlockDataWriter()
    writer.SetInput(mblock)
    writer.SetFileName('bob.vtm')
    writer.Write()

    return mblock

def fluidity_to_ugrid_p1(state, p1_check, p0_check, coordinates):
    """ Extract fields on P0 and P1 continuous meshes from fluidity state to a vtkUnstructuredGrid data object.

    Fields on P0 meshes are stored as cell data, while data on P1 continuous meshes is stored as point data."""

    dimension=    dimension=state.vector_fields['Coordinate'].dimension

    meshes = [mesh for mesh in state.meshes.values() if p1_check(mesh,
                                                                 dimension)]
    meshes_p0 = [mesh for mesh in state.meshes.values() if p0_check(mesh,
                                                                    dimension)]

    point_data = coordinates.val

    if coordinates.dimension < 3:
        point_data = numpy.concatenate((point_data,
                                        numpy.zeros((coordinates.node_count,
                                                     3-coordinates.dimension))),
                                         axis=1)

    pts = vtk.vtkPoints()
    pts.SetDataTypeToDouble()
    pts.SetNumberOfPoints(coordinates.node_count)
    data = vtk.vtkDoubleArray()
    data.DeepCopy(numpy_support.numpy_to_vtk(point_data))
    pts.SetData(data)


    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(pts)

    for k in range(coordinates.element_count):
        id_list = vtk.vtkIdList()
        shape = coordinates.ele_shape(k)

        for node in coordinates.ele_nodes(k):
            id_list.InsertNextId(node)

        ugrid.InsertNextCell(CELL_DICT[(shape.type, shape.dimension, shape.loc)],id_list)

    ugrid = fluidity_data_to_ugrid(state, meshes, ugrid)
    ugrid = fluidity_cell_data_to_ugrid(state, meshes_p0, ugrid)

    return ugrid


def is_p0(mesh, dimension):
    """ Test if mesh is a P0 mesh."""

    return (mesh.continuity == -1 and
            mesh.shape.dimension == dimension and
            mesh.shape.type == 'lagrangian' and
            mesh.shape.degree == 0)

def is_p1(mesh, dimension):
    """ Test if mesh is a P1 continous mesh."""

    return (mesh.continuity > -1 and
            mesh.shape.dimension == dimension and
            mesh.shape.type == 'lagrangian' and
            mesh.shape.degree == 1)

def is_surface_p1(mesh, dimension):
    """ Test if mesh is a P1 continous surface mesh."""

    return (mesh.continuity > -1 and
            mesh.shape.dimension == dimension-1 and
            mesh.shape.type == 'lagrangian' and
            mesh.shape.degree == 1)

def is_surface_p0(mesh, dimension):
    """ Test if mesh is a P0 mesh."""

    return (mesh.continuity == -1 and
            mesh.shape.dimension == dimension-1 and
            mesh.shape.type == 'lagrangian' and
            mesh.shape.degree == 0)

def is_p1dg(mesh, dimension):
    """ Test if mesh is a P1 discontinous mesh."""

    return (mesh.continuity == -1 and
            mesh.shape.dimension == dimension and
            mesh.shape.type == 'lagrangian' and
            mesh.shape.degree == 1)

def is_p2(mesh, dimension):
    """ Test if mesh is a P2 continuous mesh."""

    return (mesh.continuity > -1 and
            mesh.shape.dimension == dimension and
            mesh.shape.type == 'lagrangian' and
            mesh.shape.degree == 2)


def is_p2dg(mesh, dimension):
    """ Test if mesh is a P2 discontinuous mesh."""

    return (mesh.continuity == -1 and
            mesh.shape.dimension == dimension and
            mesh.shape.type == 'lagrangian' and
            mesh.shape.degree == 2)

def fluidity_to_ugrid_by_mesh(state, test):
    """ Extract fludity data on a generic mesh to a vtkUnstructuredGrid data object.

    test should be a function which returns True when passed the desired mesh type."""

    dimension=state.vector_fields['Coordinate'].dimension

    meshes = [mesh for mesh in state.meshes.values() if test(mesh, dimension)]

    if not meshes:
        return None

    coordinates = state.vector_fields['Coordinate']

    pts = vtk.vtkPoints()
    pts.Allocate(meshes[0].node_count)
    pts.SetNumberOfPoints(meshes[0].node_count)
     
    ugrid = vtk.vtkUnstructuredGrid()

    for k in range(meshes[0].element_count):
        id_list = vtk.vtkIdList()
        id_list.Allocate(meshes[0].ele_loc(k))
        id_list.SetNumberOfIds(meshes[0].ele_loc(k))
        shape = meshes[0].shape

        lpoint = numpy.zeros(3)

        for loc, (point, node) in enumerate(zip(coordinates.remap_ele(k, meshes[0]),
                                                meshes[0].ele_nodes(k))):

            lpoint[:len(point)] = point 
            pts.SetPoint(node, *lpoint)
            id_list.SetId(NUM_DICT[(shape.type, shape.dimension, shape.loc)][loc], node)

        ugrid.InsertNextCell(CELL_DICT[(shape.type, shape.dimension, shape.loc)], id_list)

    ugrid.SetPoints(pts)
    ugrid = fluidity_data_to_ugrid(state, meshes, ugrid)

    return ugrid

def fluidity_data_to_ugrid(state, meshes, ugrid):
    """ Extract fluidity data from meshes on a desired type to an existing unstructured grid object's point data."""

    for name, field in (state.scalar_fields.items()
                        +state.vector_fields.items()
                        +state.tensor_fields.items()):
        
        if field.mesh not in meshes:
            continue

        if field.val.shape[0] == 1:
            val = field.node_val(0)
            data = vtk.vtkDoubleArray()
            data.SetNumberOfComponents(numpy.prod(val.shape))
            for k in range(ugrid.GetNumberOfPoints()):
                data.InsertNextTuple(val.ravel())
        else:
            try:
                data = numpy_support.numpy_to_vtk(field.val)
            except:
                val = field.node_val(0)
                data = vtk.vtkDoubleArray()
                data.SetNumberOfComponents(numpy.prod(val.shape))
                for k in range(ugrid.GetNumberOfPoints()):
                    data.InsertNextTuple(field.node_val(k).ravel())

        data.SetName(name)
        ugrid.GetPointData().AddArray(data)

        data.SetName(name)
        ugrid.GetPointData().AddArray(data)
        
    return ugrid

def fluidity_cell_data_to_ugrid(state, meshes, ugrid):
    """ Extract fluidity data from meshes on a desired type to an existing unstructured grid object's cell data.

    Will only work for P0 meshes."""

    for name, field in (state.scalar_fields.items()
                        +state.vector_fields.items()
                        +state.tensor_fields.items()):
        
        if field.mesh not in meshes:
            continue

        if field.val.shape[0] == 1:
            val = field.node_val(0)
            data = vtk.vtkDoubleArray()
            data.SetNumberOfComponents(numpy.prod(val.shape))
            for k in range(ugrid.GetNumberOfPoints()):
                data.InsertNextTuple(val.ravel())
        else:
            data = numpy_support.numpy_to_vtk(field.val)

        data.SetName(name)
        ugrid.GetCellData().AddArray(data)
        
    return ugrid

def extract_field_from_ugrid(ugrid,field,name=None):
    """ Extract data from an unstructured grid onto a field """

    if name is None:
        name = field.name

    if is_p0(field.mesh):
        data = ugrid.GetCellData().GetArray(name)
        if data:
            ndata = numpy_support.vtk_to_numpy(data)
            assert ndata.shape == field.val.shape
            field.val[:] = ndata[:]
        else:
            print ("P0 vtk field %s not found"%name)

    else:
        data = ugrid.GetPointData().GetArray(name)
        if data:
            ndata = numpy_support.vtk_to_numpy(data)
            assert ndata.shape == field.val.shape
            field.val[:]=numpy_support.vtk_to_numpy(data)[:]
        else:
            print ("vtk field %s not found"%name)


def write_ugrid(data,basename,path='.'):
    """ Parallel safe data object writer"""

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if size == 1:

        ### serial case

        filename = basename +'.vtu'

        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(path+'/'+filename)
        writer.SetInput(data)
        writer.Write()

        return

    ## In parallel we make a directory and dump the files into it

    try:
        os.makedirs(path+'/'+basename)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    filename=basename+'.pvtu'

    writer = vtk.vtkXMLPUnstructuredGridWriter()
    writer.SetNumberOfPieces(size)
    writer.SetStartPiece(rank)
    writer.SetEndPiece(rank)

    writer.SetFileName(path+'/'+basename+'/'+filename)
    writer.SetInput(data)
    writer.Write()


    if rank == 0:
        os.rename(path+'/'+basename+'/'+filename, path+'/'+filename)

        stream=open(path+'/'+filename).read()
        stream=stream.replace('<Piece Source="'+basename, 
                              '<Piece Source="'+basename+'/'+basename)
        outfile=open(path+'/'+filename, 'w')
        outfile.write(stream)
        outfile.flush()
        outfile.close()

    return 
