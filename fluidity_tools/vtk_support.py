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

    return ugrid
        
    



