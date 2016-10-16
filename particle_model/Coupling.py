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
    if vtk.vtkVersion.GetVTKMajorVersion()<6:
        cvs.SetInput(linear_data)
    else:
        cvs.SetInputData(linear_data)
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
    if vtk.vtkVersion.GetVTKMajorVersion()<6:
        cvs.SetInput(linear_data)
    else:
        cvs.SetInputData(linear_data)
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

def barocentric_id(cell, pos):
    p0 = numpy.array(cell.GetPoints().GetPoint(0))
    p1 = numpy.array(cell.GetPoints().GetPoint(1))
    if cell.IsA('vtkLine'):
        if sum((pos-p0)**2)/sum((p1-p0)**2)>0.25:
            return cell.GetPointIds().GetId(1)
        else:
            return cell.GetPointIds().GetId(0)
    else:
        p2 = numpy.array(cell.GetPoints().GetPoint(2))
        d11 = numpy.dot(p1-p0,p1-p0)
        d12 = numpy.dot(p1-p0,p2-p0)
        d22 = numpy.dot(p2-p0,p2-p0)
        re1 =numpy.dot(pos-p0,p1-p0)
        re2 = numpy.dot(pos-p0,p2-p0)
        det = d11 * d22 - d12**2
        res = numpy.empty(3,float)
        res[0] = (d22*re1-d12*re2)/det
        res[1] = (d11*re2-d12*re1)/det
        res[2] = 1.0 - res[0] - res[1]

        return cell.GetPointIds().GetId(res.argmax())


def get_wear_rate_source(bucket, alpha, delta_t):
    """Calculate wear_rate on the boundary surface"""

    linear_data = bucket.system.boundary.bnd
    is2d = linear_data.GetCell(0).GetCellType()==vtk.VTK_LINE

    if is2d:
        dim=2
    else:
        dim=3
    wear = numpy.zeros((linear_data.GetNumberOfPoints()))
    volume = numpy.zeros(linear_data.GetNumberOfPoints())

    for _ in range(linear_data.GetNumberOfCells()):
        
        cell = linear_data.GetCell(_)
        pntIds = cell.GetPointIds()
        cv_mass = IO.get_measure(cell)/cell.GetNumberOfPoints()

        for dummy_1 in range(pntIds.GetNumberOfIds()):
            volume[pntIds.GetId(dummy_1)] += cv_mass

    for col in bucket.collisions():
        if col.time<bucket.time-bucket.delta_t:
            continue

        cell = linear_data.GetCell(col.cell)

        gid = barocentric_id(cell,col.pos)
                                                      
        wear[gid] += col.get_wear()/delta_t

    wear /= volume

            
    return wear
