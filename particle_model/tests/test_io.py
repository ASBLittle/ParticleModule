"""Unit tests for the file IO"""
from particle_model import IO

import vtk

def test_clean_unstructured_grid():
    """ Test the clean_unstructured_grid routine"""
    reader = vtk.vtkXMLUnstructuredGridReader()
    del reader
    ugrid = vtk.vtkUnstructuredGrid()
    clean_grid = IO.clean_unstructured_grid(ugrid)

    assert clean_grid


