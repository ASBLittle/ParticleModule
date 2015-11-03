""" Module to convert from the state_type to vtk. """

import vtk
from vtk.util import numpy_support
import numpy

NUM_DICT = {('lagrangian', 2, 3) : [0,1,2],
             ('lagrangian', 2, 6) : [0,3,1,5,4,2],
             ('lagrangian', 3, 4) : [0,1,2,3],
             ('lagrangian', 3, 10) : [0,]}

CELL_DICT = {('lagrangian', 2, 3) : vtk.VTK_TRIANGLE,
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

    ugrid = fluidity_to_ugrid_p1(state)
    mblock.SetBlock(0,ugrid)

    dummy = 1

    for mesh_test in (is_p1dg, is_p2, is_p2dg):
        ugrid = fluidity_to_ugrid_by_mesh(state, mesh_test)
        if ugrid:
            mblock.SetBlock(dummy, ugrid)
            dummy += 1

    return mblock

def fluidity_to_ugrid_p1(state):
    """ Extract fields on P0 and P1 continuous meshes from fluidity state to a vtkUnstructured Grid data object."""

    meshes = [mesh for mesh in state.meshes.values() if is_p1(mesh)]
    meshes_p0 = [mesh for mesh in state.meshes.values() if is_p0(mesh)]

    coordinates = state.vector_fields['Coordinate']

    point_data = coordinates.val

    if coordinates.dimension < 3:
        point_data = numpy.concatenate((point_data,
                                        numpy.zeros((coordinates.node_count,
                                                     3-coordinates.dimension))),
                                         axis=1)

    pts = vtk.vtkPoints()
    pts.SetData(numpy_support.numpy_to_vtk(point_data))
    pts.SetNumberOfPoints(coordinates.node_count)

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


def is_p0(mesh):
    """ Test if mesh is a P0 mesh."""

    return (mesh.continuity > -1 and
            mesh.shape.type == 'lagrangian' and
            mesh.shape.degree == 0)

def is_p1(mesh):
    """ Test if mesh is a P1 continous mesh."""

    return (mesh.continuity > -1 and
            mesh.shape.type == 'lagrangian' and
            mesh.shape.degree == 1)

def is_p1dg(mesh):
    """ Test if mesh is a P1 discontinous mesh."""

    return (mesh.continuity == -1 and
            mesh.shape.type == 'lagrangian' and
            mesh.shape.degree == 1)

def is_p2(mesh):
    """ Test if mesh is a P2 continuous mesh."""

    return (mesh.continuity > -1 and
            mesh.shape.type == 'lagrangian' and
            mesh.shape.degree == 2)

def is_p2dg(mesh):
    """ Test if mesh is a P2 discontinuous mesh."""

    return (mesh.continuity == -1 and
            mesh.shape.type == 'lagrangian' and
            mesh.shape.degree == 2)

def fluidity_to_ugrid_by_mesh(state, test):
    """ Extract fludity data on a generic mesh to a vtkUnstructuredGrid data object.

    test should be a function which returns True when passed the desired mesh type."""

    meshes = [mesh for mesh in state.meshes.values() if test(mesh)]

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
            data = numpy_support.numpy_to_vtk(field.val)

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


        
    



