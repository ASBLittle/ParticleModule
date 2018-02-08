import vtk
import vtk_extras
import numpy

reader = vtk.vtkXMLGenericDataObjectReader()
reader.SetFileName('../examples/lagrangian/gyre40x40_0.vtu')
reader.Update()
ugrid = reader.GetOutput()

locator=vtk.vtkCellLocator()
locator.SetDataSet(ugrid)
pick=vtk_extras.Picker(locator)
pick.locator = locator
print pick.locator
print pick((1.0,0.5,0.0), ugrid.GetPointData().GetArray("Velocity"))
