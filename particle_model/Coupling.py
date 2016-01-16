""" Module dealing with coupling in fluidity"""

import particle_model.vtkParticlesPython as vtp
from particle_model import IO

import vtk
import numpy

ARGV = [0.0, 0.0, 0.0]
ARGI = vtk.mutable(0)
ARGR = vtk.mutable(0.0)

def get_cv_fraction(bucket, data):
    """Calculate the particle volume fraction using control volumes"""

    linear_data = IO.get_linear_block(data) 
    is2d = linear_data.GetCell(0).GetCellType()==vtk.VTK_TRIANGLE

    cvs = vtp.vtkShowCVs()
    cvs.SetContinuity(-1)
    cvs.SetInput(linear_data)
    cvs.Update()
    cv_data = cvs.GetOutput()

    locator = vtk.vtkCellLocator()
    locator.SetDataSet(cv_data)
    locator.BuildLocator()

    output = numpy.zeros(linear_data.GetNumberOfPoints())
    volume = numpy.zeros(linear_data.GetNumberOfPoints())

    for _ in range(linear_data.GetNumberOfCells()):
        
        cell = linear_data.GetCell(_)
        pntIds = cell.GetPointIds()
        cv_mass = IO.get_measure(cell)/cell.GetNumberOfPoints()

        for dummy_1 in range(pntIds.GetNumberOfIds()):
            volume[pntIds.GetId(dummy_1)] += cv_mass

    for par in bucket.particles:
        index = locator.FindCell(par.pos)
        if index<0:
            continue
        ele, l_id = divmod(index,linear_data.GetCell(0).GetNumberOfPoints())

        gid = linear_data.GetCell(ele).GetPointId(l_id)

                                                      
        if is2d:
            output[gid] += par.parameters.get_area()
        else:
            output[gid] += par.parameters.get_volume()
    


    return output/volume


def get_solid_velocity(bucket, data, volfrac):
    """Calculate the particle volume fraction using control volumes"""

    linear_data = IO.get_linear_block(data) 
    is2d = linear_data.GetCell(0).GetCellType()==vtk.VTK_TRIANGLE

    cvs = vtp.vtkShowCVs()
    cvs.SetContinuity(-1)
    cvs.SetInput(linear_data)
    cvs.Update()
    cv_data = cvs.GetOutput()

    locator = vtk.vtkCellLocator()
    locator.SetDataSet(cv_data)
    locator.BuildLocator()


    if is2d:
        dim=2
    else:
        dim=3
    output = numpy.zeros((linear_data.GetNumberOfPoints(),dim))
    volume = numpy.zeros(linear_data.GetNumberOfPoints())


    for par in bucket.particles:
        index = locator.FindCell(par.pos)
        if index<0:
            continue
        ele, l_id = divmod(index,linear_data.GetCell(0).GetNumberOfPoints())

        gid = linear_data.GetCell(ele).GetPointId(l_id)

                                                      
        if is2d:
            volume[gid] += par.parameters.get_area()
            output[gid,:] += par.parameters.get_area()*par.vel[:dim]
        else:
            volume[gid] += par.parameters.get_volume()
            output[gid,:] += par.parameters.get_volume()*par.vel[:dim]

    for _ in range(linear_data.GetNumberOfPoints()):
        if volume[_]>0.0:
            output[_,:] = output[_,:]/volume[_]

            
    return output
