import vtk
from particle_model import vtk_extras

reader = vtk.vtkXMLGenericDataObjectReader()

reader.SetFileName('particle_model/tests/data/gyre_0.vtu')
reader.Update()
ugrid = reader.GetOutput()

def test_BoundingSurface():
    """Test the vtk_extras.BoundingSurface function."""
    pass

def test_FindCell():
    """Test the vtk_extras.FindCell function."""
    locator = vtk.vtkCellLocator()
    locator.SetDataSet(ugrid)

def test_EvaluateField():
    """Test the vtk_extras.EvaluateField function."""
    pass

def test_vInterpolate():
    """Test the vtk_extras.EvaluateField function."""
    pass

def test_Picker():
    """Test the vtk_extras.Picker class"""
    locator = vtk.vtkCellLocator()
    locator.SetDataSet(ugrid)
    picker = vtk_extras.Picker()
    picker.name = "Velocity"
    picker.grid = ugrid
    picker.locator = locator

#    out = picker((0.5,0.5,0.0))
