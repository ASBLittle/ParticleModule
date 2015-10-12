""" Module containing input-output routines between the particle model and
the file system. Mostly vtk."""

import Collision

import vtk
import pylab as p
import numpy
import os
import os.path


types_3d = [vtk.VTK_TETRA, vtk.VTK_QUADRATIC_TETRA]
types_2d = [vtk.VTK_TRIANGLE, vtk.VTK_QUADRATIC_TRIANGLE]
types_1d = [vtk.VTK_LINE]

class boundaryData(object):
    """ Class storing the boundary data for the problem"""
    def __init__(self, filename):
        """Class containing the information about the boundary of the domain.

        Args:
            filename (str): Name of the file containing the
            vtkUnstructuredGrid denoting the boundary of the domain."""

        self.geom_filter = vtk.vtkGeometryFilter()
        self.reader = vtk.vtkXMLUnstructuredGridReader()
        self.bnd = self.reader.GetOutput()
        self.bndl = vtk.vtkCellLocator()

        self.update_boundary_file(filename)

    def update_boundary_file(self, filename):
        """ Update the boundary data from the file."""
        if not os.path.isfile(filename):
            print os.getcwd()
            raise OSError

        self.reader.SetFileName(filename)
        self.reader.Update()
        self.bnd.Update()


        self.geom_filter.SetInput(self.bnd)
        self.geom_filter.Update()

        self.bndl.SetDataSet(self.geom_filter.GetOutput())
        self.bndl.BuildLocator()

    def rebuild_locator(self):
        """ Rebuild the locator information"""
        self.bndl.BuildLocator()



def clean_unstructured_grid(ugrid):
    """Collapse a vtu produced from a discontinuous grid back down to the continuous space.

    Args:
    ugrid (vtkUnstructuredGrid): the input discontinuous grid

    Results
    out_grid (vtkUnstructuredGrid): A continuous grid"""

    merge_points = vtk.vtkMergePoints()
    out_grid = vtk.vtkUnstructuredGrid()

    for i in range(ugrid.GetNumberOfPoints()):
        merge_points.InsertUniquePoint(ugrid.GetPoints().GetPoint(i))

    merge_points.BuildLocator()

    pts = vtk.vtkPoints()
    pts.DeepCopy(merge_points.GetPoints())
    out_grid.SetPoints(pts)

    for i in range(ugrid.GetNumberOfCells()):
        cell = ugrid.GetCell(i)
        cell_ids = cell.GetPointIds()

        for j in range(cell.GetNumberOfPoints()):

            original_point = cell.GetPoints().GetPoint(j)
            cell_ids.SetId(j,
                           merge_points.FindClosestInsertedPoint(original_point))

        out_grid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())


    out_grid.GetCellData().DeepCopy(ugrid.GetCellData())

    return out_grid

def extract_boundary(ugrid):
    """Extract the boundary elements from an unstructured grid, provided it already contains them.

    Args:

    ugrid (vtkUnstructuredGrid): The grid with which to work.

    Results:

    out_grid (vtkUnstructuredGrid): Grid containing the boundary of ugrid"""

    out_grid = vtk.vtkUnstructuredGrid()
    pts = vtk.vtkPoints()
    pts.DeepCopy(ugrid.GetPoints())
    out_grid.SetPoints(pts)
    out_grid.GetCellData().CopyStructure(ugrid.GetCellData())

    celltypes = vtk.vtkCellTypes()

    ugrid.GetCellTypes(celltypes)

    if any([celltypes.IsType(ct) for ct in types_3d]):
        dim = 3
    elif any([celltypes.IsType(ct) for ct in types_2d]):
        dim = 2
    elif any([celltypes.IsType(ct) for ct in types_1d]):
        dim = 1
    else:
        dim = 0

    print dim

    NC = ugrid.GetNumberOfCells()
    NCDA = ugrid.GetCellData().GetNumberOfArrays()

    for i in range(NCDA):
        out_grid.GetCellData().GetArray(i).SetName(ugrid.GetCellData().GetArray(i).GetName())

    cell_data = ugrid.GetCellData()
    for i in range(NC):
        cell = ugrid.GetCell(i)
        if dim > cell.GetCellDimension():
            out_grid.InsertNextCell(cell.GetCellType(),
                                    cell.GetPointIds())
            for j in range(NCDA):
                out_data = out_grid.GetCellData().GetArray(j)
                out_data.InsertNextTuple(cell_data.GetArray(j).GetTuple(i))

    return out_grid

def plot_boundary(ugrid, **kwargs):

    """Plot a boundary using matplotlib

    Args:
        ugrid (vtkUnstructuredGrid): The boundary mesh to plot

    Other arguments are passed on to the matplotlib plot command"""

    for i in range(ugrid.GetNumberOfCells()):
        cell = ugrid.GetCell(i)

        pos_x = []
        pos_y = []

        for j in range(cell.GetNumberOfPoints()):
            pnt = cell.GetPoints().GetPoint(j)
            pos_x.append(pnt[0])
            pos_y.append(pnt[1])

        p.plot(pos_x, pos_y, 'k', **kwargs)

def get_ascii_data(filename='data.dat'):

    """Read the ascii output file and output as numpy data arrays

    Args:
        filename (str) The file to interogate

    Results:
        t (ndarray): Time
        x (ndarray): x coordinate of particle position
        y (ndarray): y coordinate of particle position
        z (ndarray): z coordinate of particle position
        u (ndarray): u coordinate of particle velocity
        v (ndarray): v coordinate of particle velocity
        w (ndarray): w coordinate of particle velocity"""

    infile = open(filename, 'r')

    time = []
    x = []
    y = []
    z = []
    u = []
    v = []
    w = []

    for line in infile.readlines():
        data = [float(M) for M in line.split()]

        time.append(data[0])
        n = len(data[1:])/6
        x.append(data[1::3][:n])
        u.append(data[1::3][n:])
        y.append(data[2::3][:n])
        v.append(data[2::3][:n])
        z.append(data[3::3][:n])
        w.append(data[3::3][n:])

    x = numpy.array(x)
    y = numpy.array(y)
    z = numpy.array(z)
    u = numpy.array(u)
    v = numpy.array(v)
    w = numpy.array(w)

    infile.close()

    return time, x, y, z, u, v, w


def ascii_to_polydata_time_series(filename,basename):

    """Convert ascii file to a series of vtkPolyData (.vtp) files. 

    Each file contains one time level of the data, and are numbered sequentially. Within each file, each dataset is written to seperate pixel.

    Args:
        filename (str): Filename/path of the ascii file containing the data.
        basename (str): String used in the construction of the file series. The formula is of the form basename_0.vtp, basename_1.vtp,..."""

    time, pos_x, pos_y, pos_z, vel_u, vel_v, vel_w = get_ascii_data(filename)

    for i, full_data in enumerate(zip(pos_x, pos_y, pos_z,
                                               vel_u, vel_v, vel_w)):
        poly_data = vtk.vtkPolyData()
        pnts = vtk.vtkPoints()
        pnts.Allocate(0)
        poly_data.SetPoints(pnts)
        poly_data.Allocate(pos_x.shape[0])

        outtime = vtk.vtkDoubleArray()
        outtime.Allocate(pos_x.shape[0])
        outtime.SetName('Time')

        velocity = vtk.vtkDoubleArray()
        velocity.SetNumberOfComponents(3)
        velocity.Allocate(pos_x.shape[0])
        velocity.SetName('Particle Velocity')

        for k, data in enumerate(full_data):
            pixel = vtk.vtkPixel()
            pixel.GetPointIds().InsertId(k,
                                         poly_data.GetPoints().InsertNextPoint(data[0],
                                                                               data[1],
                                                                               data[2]))
            outtime.InsertNextValue(time[i])
            velocity.InsertNextTuple3(data[3], data[4], data[5])

        poly_data.GetPointData().AddArray(outtime)
        poly_data.GetPointData().AddArray(velocity)
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName("%s_%d.vtp"%(basename, i))
        writer.SetInput(poly_data)
        writer.Write()

def ascii_to_polydata(filename, outfile):
    """Convert ascii file to a single vtkPolyData (.vtp) files.

    Each particle is written to seperate cell.

    Args:
        filename (str): Filename/path of the ascii file containing the data.
        outfile (str):  Filename of the output PolyDataFile. The extension .vtp
        is NOT added automatically."""

    poly_data = vtk.vtkPolyData()
    pnts = vtk.vtkPoints()
    pnts.Allocate(0)
    poly_data.SetPoints(pnts)
    full_data = get_ascii_data(filename)
    time = full_data[0]
    pos_x = full_data[1]
    pos_y = full_data[2]
    pos_z = full_data[3]
    poly_data.Allocate(pos_x.shape[1])

    outtime = vtk.vtkDoubleArray()
    outtime.SetName('Time')

    for X, Y, Z in zip(pos_x.T, pos_y.T, pos_z.T):
        line = vtk.vtkLine()
        for k, data in enumerate(zip(time, X, Y, Z)):
            line.GetPointIds().InsertId(k,
                                        poly_data.GetPoints().InsertNextPoint(data[1],
                                                                              data[2],
                                                                              data[3]))
            outtime.InsertNextValue(data[0])
        poly_data.InsertNextCell(line.GetCellType(), line.GetPointIds())

    poly_data.GetPointData().AddArray(outtime)
    
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(outfile)
    writer.SetInput(poly_data)

    writer.Write()

def collision_list_to_polydata(col_list, outfile,
                               model=Collision.MclauryMassCoeff, **kwargs):
    """Convert collision data to a single vtkPolyData (.vtp) files. 

    Each particle is written to seperate cell.

    Args:
        filename (str): Filename/path of the ascii file containing the data.
        outfile (str):  Filename of the output PolyDataFile. The extension .vtp
        is NOT added automatically."""

    poly_data = vtk.vtkPolyData()
    pnts = vtk.vtkPoints()
    pnts.Allocate(0)
    poly_data.SetPoints(pnts)
    poly_data.Allocate(len(col_list))

    time = vtk.vtkDoubleArray()
    time.SetName('Time')
    wear = vtk.vtkDoubleArray()
    wear.SetName('Wear')
        
    for col in col_list:
        pixel = vtk.vtkPixel()
        pixel.GetPointIds().InsertId(0, 
                                     poly_data.GetPoints().InsertNextPoint(col.x[0],
                                                                           col.x[1],
                                                                           col.x[2]))
        time.InsertNextValue(col.time)
        wear.InsertNextValue(model(col, **kwargs))
        poly_data.InsertNextCell(pixel.GetCellType(), pixel.GetPointIds())

    poly_data.GetPointData().AddArray(time)
    poly_data.GetPointData().AddArray(wear)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(outfile)
    writer.SetInput(poly_data)

    writer.Write()
