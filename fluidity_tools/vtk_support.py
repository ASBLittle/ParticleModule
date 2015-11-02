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

    pts = vtk.vtkPoints()
    pts.Allocate(coordinates.node_count)

    for point in coordinates.val:
        lpoint = numpy.zeros(3)
        lpoint[:len(point)] = point 
        pts.InsertNextPoint(*lpoint)

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

        for loc, (point, node) in enumerate(zip(coordinates.remap_ele(k, meshes[0]),
                                                meshes[0].ele_nodes(k))):

            lpoint = numpy.zeros(3)
            lpoint[:len(point)] = point 
            pts.SetPoint(node, *lpoint)
            id_list.SetId(NUM_DICT[(shape.type, shape.dimension, shape.loc)][loc], node)

        ugrid.InsertNextCell(CELL_DICT[(shape.type, shape.dimension, shape.loc)], id_list)

    ugrid.SetPoints(pts)
    ugrid = fluidity_data_to_ugrid(state, meshes, ugrid)

    return ugrid

def fluidity_data_to_ugrid(state, meshes, ugrid):
    """ Extract fluidity data from meshes on a desired type to an existing unstructured grid object's point data."""

    for name, scalar in state.scalar_fields.items():
        
        if scalar.mesh not in meshes:
            continue

        if scalar.node_count == 1:
            data = vtk.vtkDoubleArray()
            data.Allocate(ugrid.GetNumberOfPoints())
            for k in range(ugrid.GetNumberOfPoints()):
                data.InsertNextValue(scalar[0])
        else:
            data = numpy_support.numpy_to_vtk(scalar.val)

        data.SetName(name)
        ugrid.GetPointData().AddArray(data)

    for name, vector in state.vector_fields.items():
        
        if vector.mesh not in meshes:
            continue

        data = vtk.vtkDoubleArray()
        data.SetNumberOfComponents(3)

        if vector.node_count == 1:
            for k in range(ugrid.GetNumberOfPoints()):
                lval=numpy.zeros(3)
                lval[:vector.val.shape[1]] = vector.node_val(0)
                data.InsertNextTuple3(*(lval))
        else:
            for k in range(ugrid.GetNumberOfPoints()):
                lval=numpy.zeros(3)
                lval[:vector.val.shape[1]] = vector.node_val(k)
                data.InsertNextTuple3(*(lval))

        data.SetName(name)
        ugrid.GetPointData().AddArray(data)


    for name, tensor in state.tensor_fields.items():
        
        if tensor.mesh not in meshes:
            continue

        data = vtk.vtkDoubleArray()
        data.SetNumberOfComponents(9)

        if tensor.val.shape[0] == 1:
            for k in range(ugrid.GetNumberOfPoints()):
                lval=numpy.zeros((3, 3))
                lval[:tensor.val.shape[1],:tensor.val.shape[2]] = tensor.node_val(0)
                data.InsertNextTuple9(*(lval.ravel()))
        else:
            for k in range(ugrid.GetNumberOfPoints()):
                lval=numpy.zeros((3, 3))
                lval[:tensor.val.shape[1],:tensor.val.shape[2]] = tensor.node_val(k)
                data.InsertNextTuple9(*(lval.ravel()))

        data.SetName(name)
        ugrid.GetPointData().AddArray(data)
        
    return ugrid

def fluidity_cell_data_to_ugrid(state, meshes, ugrid):
    """ Extract fluidity data from meshes on a desired type to an existing unstructured grid object's cell data.

    Will only work for P0 meshes."""

    for name, scalar in state.scalar_fields.items():
        
        if scalar.mesh not in meshes:
            continue

        if scalar.node_count == 1:
            data = vtk.vtkDoubleArray()
            data.Allocate(ugrid.GetNumberOfCells())
            for k in range(ugrid.GetNumberOfCells()):
                data.InsertNextValue(scalar[0])
        else:
            data = numpy_support.numpy_to_vtk(scalar.val)

        data.SetName(name)
        ugrid.GetCellData().AddArray(data)

    for name, vector in state.vector_fields.items():
        
        if vector.mesh not in meshes:
            continue

        data = vtk.vtkDoubleArray()
        data.SetNumberOfComponents(3)

        if vector.node_count == 1:
            for k in range(ugrid.GetNumberOfCells()):
                lval=numpy.zeros(3)
                lval[:vector.val.shape[1]] = vector.node_val(0)
                data.InsertNextTuple3(*(lval))
        else:
            for k in range(ugrid.GetNumberOfCells()):
                lval=numpy.zeros(3)
                lval[:vector.val.shape[1]] = vector.node_val(k)
                data.InsertNextTuple3(*(lval))

        data.SetName(name)
        ugrid.GetCellData().AddArray(data)


    for name, tensor in state.tensor_fields.items():
        
        if tensor.mesh not in meshes:
            continue

        data = vtk.vtkDoubleArray()
        data.SetNumberOfComponents(9)

        if tensor.val.shape[0] == 1:
            for k in range(ugrid.GetNumberOfCells()):
                lval=numpy.zeros((3, 3))
                lval[:tensor.val.shape[1],:tensor.val.shape[2]] = tensor.node_val(0)
                data.InsertNextTuple9(*(lval.ravel()))
        else:
            for k in range(ugrid.GetNumberOfCells()):
                lval=numpy.zeros((3, 3))
                lval[:tensor.val.shape[1],:tensor.val.shape[2]] = tensor.node_val(k)
                data.InsertNextTuple9(*(lval.ravel()))

        data.SetName(name)
        ugrid.GetCellData().AddArray(data)
        
    return ugrid


        
    


