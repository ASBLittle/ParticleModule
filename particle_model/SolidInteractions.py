"""Module containing routines to implement solid-solid interactions via granular
 temperature."""

import numpy

from particle_model import vtkParticlesPython
from particle_model.Debug import logger

import vtk
from vtk.util import numpy_support

def radial_distribution_function(alpha):

    """Calculate radial distribution function."""

    ALPHA_MAX = 0.59999

    ALPHA0 =0.6

    numpy.where(alpha > ALPHA_MAX, ALPHA_MAX, alpha)

    return 1.0/(1.0-(alpha/ALPHA0)**(1.0/3.0))

def rdf_deriv(alpha):

    """Derivative of radial distribution function."""

    ALPHA_MAX = 0.59999

    ALPHA0 = 0.6

    numpy.where(alpha > ALPHA_MAX, ALPHA_MAX, alpha)

    return -1.0/(3.0*ALPHA0)*(alpha/ALPHA0)**(-2.0/3.0)/(1.0-(alpha/ALPHA0)**(1.0/3.0))**2

def distance2(pnt1, pnt2):
    """Get square of distance between two points."""
    return vtk.vtkMath().Distance2BetweenPoints(pnt1, pnt2)

def calculate_averaged_properties_cpp(poly_data):
    """Calculate volume averaged values using C++."""

    gt_filter = vtkParticlesPython.vtkGranularTemperature()
    if vtk.vtkVersion.GetVTKMajorVersion() < 6:
        gt_filter.SetInput(poly_data)
    else:
        gt_filter.SetInputData(poly_data)
    gt_filter.Update()
    poly_data.DeepCopy(gt_filter.GetOutput())

    return numpy_support.vtk_to_numpy(poly_data.GetPointData().GetVectors("SolidPressureGradient"))

def calculate_averaged_properties(poly_data, bucket):

    """ Calculate a volume fraction estimate at the level of the particles."""

    locator = vtk.vtkPointLocator()
    locator.SetDataSet(poly_data)
    locator.BuildLocator()

    LENGTH = 0.03
    MODIFIER = 3e3

    volume = numpy.zeros(poly_data.GetNumberOfPoints())
    temperature = numpy.zeros(poly_data.GetNumberOfPoints())
    solid_pressure = numpy.zeros(poly_data.GetNumberOfPoints())
    velocity = numpy.zeros((poly_data.GetNumberOfPoints(), 3))
    solid_pressure_gradient = numpy.zeros((poly_data.GetNumberOfPoints(), 3))

    for particle in bucket:
        point_list = vtk.vtkIdList()
        locator.FindPointsWithinRadius(LENGTH, particle.pos, point_list)

        beta = 1.0/6.0*numpy.pi*particle.parameters.diameter**3

        for _ in range(point_list.GetNumberOfIds()):
            point_index = point_list.GetId(_)

            particle2 = bucket.particles[point_index]

            rad2 = distance2(particle2.pos, particle.pos)
            rad2 /= LENGTH**2

            gamma = beta*numpy.exp(-rad2)*MODIFIER

            volume[point_index] += gamma

            velocity[point_index, :] += particle.vel*gamma

    volume /= 0.5*LENGTH**2*(1.0-numpy.exp(-1.0**2))
    velocity /= 0.5*LENGTH**2*(1.0-numpy.exp(-1.0**2))

    for i in range(3):
        velocity[:, i] /= volume

    for k, particle in enumerate(bucket):
        point_list = vtk.vtkIdList()
        locator.FindPointsWithinRadius(LENGTH, particle.pos, point_list)

        beta = 1.0/6.0*numpy.pi*particle.parameters.diameter**3

        for _ in range(point_list.GetNumberOfIds()):
            point_index = point_list.GetId(_)

            rad2 = distance2(poly_data.GetPoints().GetPoint(point_index), particle.pos)
            rad2 /= LENGTH**2

            gamma = beta*numpy.exp(-rad2)*MODIFIER

            c = distance2(particle.vel, velocity[k, :])

            temperature[point_index] += c*gamma


    for particle in bucket:
        point_list = vtk.vtkIdList()
        locator.FindPointsWithinRadius(LENGTH, particle.pos, point_list)

        beta = 1.0/6.0*numpy.pi*particle.parameters.diameter**3

        for _ in range(point_list.GetNumberOfIds()):
            point_index = point_list.GetId(_)

            rad2 = distance2(poly_data.GetPoints().GetPoint(point_index), particle.pos)
            rad2 /= LENGTH **2

            gamma = beta*numpy.exp(-rad2)*MODIFIER

            c = distance2(particle.vel, velocity[point_index, :])

            val = (bucket.particles[point_index].pos-particle.pos)/LENGTH**2

            spg = ((radial_distribution_function(volume[point_index])
                    +volume[point_index]*rdf_deriv(volume[point_index]))*temperature[point_index]
                   +c*volume[point_index]*radial_distribution_function(volume[point_index]))

            solid_pressure_gradient[point_index, :] += (val*spg*gamma)

    for _ in range(poly_data.GetNumberOfPoints()):

        solid_pressure[_] = (bucket.particles[0].parameters.rho*volume[_]
                             *radial_distribution_function(volume[_])*temperature[_])

    data = [vtk.vtkDoubleArray()]
    data[0].SetName('SolidVolumeFraction')
    data.append(vtk.vtkDoubleArray())
    data[1].SetName('SolidVolumeVelocity')
    data[1].SetNumberOfComponents(3)
    data.append(vtk.vtkDoubleArray())
    data[2].SetName('GranularTemperature')
    data.append(vtk.vtkDoubleArray())
    data[3].SetName('SolidPressure')
    data.append(vtk.vtkDoubleArray())
    data[4].SetName('SolidPressureGradient')
    data[4].SetNumberOfComponents(3)

    for _ in range(poly_data.GetNumberOfPoints()):
        data[0].InsertNextValue(volume[_])
        data[1].InsertNextTuple3(*velocity[_])
        data[2].InsertNextValue(temperature[_])
        data[3].InsertNextValue(solid_pressure[_])
        data[4].InsertNextTuple3(*solid_pressure_gradient[_])

    for _ in data:
        poly_data.GetPointData().AddArray(_)

    return data[4]

def get_measure(cell):
    """Get the measure of a VTK cell."""
    if cell.GetCellType() == vtk.VTK_LINE:
        return numpy.sqrt(cell.GetLength2())
    elif cell.GetCellType() == vtk.VTK_TRIANGLE:
        return cell.ComputeArea()
    elif cell.GetCellType() == vtk.VTK_TETRA:
        pts = [cell.GetPoints().GetPoint(i) for i in range(4)]
        return cell.ComputeVolume(*pts)
    return None

def point_average(model, bucket):
    """ Calculate a volume fraction estimate at the level of the grid."""

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.DeepCopy(model)

    locator = vtk.vtkPointLocator()
    locator.SetDataSet(ugrid)
    locator.BuildLocator()

    LENGTH = 0.05

    volfrac = numpy.zeros(ugrid.GetNumberOfPoints())
    volume = numpy.zeros(ugrid.GetNumberOfPoints())
    cell_volume = numpy.zeros(ugrid.GetNumberOfPoints())
    temperature = numpy.zeros(ugrid.GetNumberOfPoints())
    solid_pressure = numpy.zeros(ugrid.GetNumberOfPoints())
    velocity = numpy.zeros((ugrid.GetNumberOfPoints(), 3))

    for _ in range(ugrid.GetNumberOfCells()):
        cell = ugrid.GetCell(_)

        loc_vol = get_measure(cell)/cell.GetNumberOfPoints()

        for i in range(cell.GetNumberOfPoints()):
            logger.info(cell.GetPointIds().GetId(i))
            cell_volume[cell.GetPointIds().GetId(i)] += loc_vol

    for particle in bucket:
        point_list = vtk.vtkIdList()
        locator.FindPointsWithinRadius(LENGTH, particle.pos, point_list)

        for _ in range(point_list.GetNumberOfIds()):
            point_index = point_list.GetId(_)

            rad2 = 0.0*distance2(ugrid.GetPoints().GetPoint(point_index), particle.pos)
            rad2 /= LENGTH**2

            gamma = particle.volume*numpy.exp(-rad2)

            volume[point_index] += gamma
            velocity[point_index, :] += particle.vel*gamma

    for _ in range(ugrid.GetNumberOfPoints()):
        if volume[_] > 1.0e-12:
            velocity[_, :] /= volume[_]

    volfrac = volume/cell_volume

    for particle in bucket:
        point_list = vtk.vtkIdList()
        locator.FindPointsWithinRadius(LENGTH, particle.pos, point_list)

        for _ in range(point_list.GetNumberOfIds()):
            point_index = point_list.GetId(_)

            rad2 = distance2(ugrid.GetPoints().GetPoint(point_index), particle.pos)
            rad2 /= LENGTH**2

            gamma = particle.volume*numpy.exp(-rad2)

            c = distance2(particle.vel, velocity[point_index, :])

            temperature[point_index] += c*gamma



    for _ in range(ugrid.GetNumberOfPoints()):
        if volume[_] > 1.0e-12:
            temperature[_] /= volume[_]

    solid_pressure = (bucket.particles[0].parameters.rho*volfrac
                      *radial_distribution_function(volfrac)*temperature)

    data = [vtk.vtkDoubleArray()]
    data[0].SetName('SolidVolumeFraction')
    data.append(vtk.vtkDoubleArray())
    data[1].SetName('SolidVolumeVelocity')
    data[1].SetNumberOfComponents(3)
    data.append(vtk.vtkDoubleArray())
    data[2].SetName('GranularTemperature')
    data.append(vtk.vtkDoubleArray())
    data[3].SetName('SolidPressure')

    for _ in range(ugrid.GetNumberOfPoints()):
        data[0].InsertNextValue(cell_volume[_])
        data[1].InsertNextTuple3(*(velocity[_]))
        data[2].InsertNextValue(temperature[_])
        data[3].InsertNextValue(solid_pressure[_])

    pdata = vtk.vtkDoubleArray()
    pdata.SetName('Time')

    for _ in range(ugrid.GetNumberOfPoints()):
        pdata.InsertNextValue(bucket.time)

    for _ in data:
        ugrid.GetPointData().AddArray(_)
    ugrid.GetPointData().AddArray(pdata)

    return ugrid

def cell_average(model, bucket):
    """ Calculate a volume fraction estimate at the level of the grid."""

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.DeepCopy(model)

    locator = vtk.vtkCellLocator()
    locator.SetDataSet(ugrid)
    locator.BuildLocator()

    volfrac = numpy.zeros(ugrid.GetNumberOfCells())
    volume = numpy.zeros(ugrid.GetNumberOfCells())
    temperature = numpy.zeros(ugrid.GetNumberOfCells())
    velocity = numpy.zeros((ugrid.GetNumberOfCells(), 3))

    for particle in bucket:
        cell_id = locator.FindCell(particle.pos)
        volume[cell_id] += particle.volume
        velocity[cell_id, :] += particle.volume*particle.vel

    for _ in range(ugrid.GetNumberOfCells()):
        if volume[_] > 1.0e-12:
            velocity[_, :] /= volume[_]
        volfrac[_] = volume[_] / get_measure(ugrid.GetCell(_))

    for particle in bucket:
        cell_id = locator.FindCell(particle.pos)
        temperature[cell_id] += particle.volume*distance2(particle.vel, velocity[cell_id, :])

    for _ in range(ugrid.GetNumberOfCells()):
        if volume[_] > 1.0e-12:
            temperature[_] /= volume[_]

    data = [vtk.vtkDoubleArray()]
    data[0].SetName('SolidVolumeFraction')
    data.append(vtk.vtkDoubleArray())
    data[1].SetName('SolidVolumeVelocity')
    data[1].SetNumberOfComponents(3)
    data.append(vtk.vtkDoubleArray())
    data[2].SetName('GranularTemperature')
#    data.append(vtk.vtkDoubleArray())
#    data[3].SetName('SolidPressure')

    for _ in range(ugrid.GetNumberOfCells()):
        data[0].InsertNextValue(volume[_])
        data[1].InsertNextTuple3(*(velocity[_]))
        data[2].InsertNextValue(temperature[_])
#        data[3].InsertNextValue(solid_pressure[_])

    pdata = vtk.vtkDoubleArray()
    pdata.SetName('Time')

    for _ in range(ugrid.GetNumberOfPoints()):
        pdata.InsertNextValue(bucket.time)

    for _ in data:
        ugrid.GetCellData().AddArray(_)
    ugrid.GetPointData().AddArray(pdata)

    return ugrid
