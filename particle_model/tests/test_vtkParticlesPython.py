import particle_model.vtkParticlesPython as vtp
import vtk


def test_showCVs():
    rdr = vtk.vtkXMLUnstructuredGridReader()
    rdr.SetFileName('particle_model/tests/data/Circle_0.vtu')
    rdr.Update()
    ugrid = rdr.GetOutput()

    cvs = vtp.vtkShowCVs()
#    cvs.SetContinuity(1)
    if vtk.VTK_MAJOR_VERSION<6:
        cvs.SetInput(ugrid)
    else:
        cvs.SetInputData(ugrid)
    cvs.Update()

#    out = cvs.GetOutput()

#    print out.GetNumberOfCells(), ugrid.GetNumberOfPoints()
    
#    assert out.GetNumberOfCells()==ugrid.GetNumberOfPoints()
