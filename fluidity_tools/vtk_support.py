""" Module to convert from the state_type to vtk. """

import vtk
import numpy

CELL_DICT = {('lagrangian', 2, 3) : vtk.VTK_TRIANGLE,
             ('lagrangian', 3, 4) : vtk.VTK_TETRA}

def fluidity_to_ugrid(state):

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

    for name, scalar in state.scalar_fields.items():

        data = vtk.vtkDoubleArray()
        data.SetName(name)
        data.Allocate(ugrid.GetNumberOfPoints())

        if scalar.node_count == 1:
            for k in range(ugrid.GetNumberOfPoints()):
                data.InsertNextValue(scalar[0])
        else:
            for k in range(ugrid.GetNumberOfPoints()):
                data.InsertNextValue(scalar[k])

        print name, ugrid.GetNumberOfPoints(), data.GetNumberOfTuples()
        ugrid.GetPointData().AddArray(data)

    for name, vector in state.vector_fields.items():
        
        data = vtk.vtkDoubleArray()
        data.SetName(name)
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

        print name, ugrid.GetNumberOfPoints(), data.GetNumberOfTuples()
        ugrid.GetPointData().AddArray(data)



    return ugrid


        
    



