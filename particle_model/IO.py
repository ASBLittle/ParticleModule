""" Module containing input-output routines between the particle model and
the file system. Mostly vtk."""

from particle_model import Collision

import vtk
import numpy
import os
import os.path
from scipy.interpolate import griddata

TYPES_3D = [vtk.VTK_TETRA, vtk.VTK_QUADRATIC_TETRA]
TYPES_2D = [vtk.VTK_TRIANGLE, vtk.VTK_QUADRATIC_TRIANGLE]
TYPES_1D = [vtk.VTK_LINE]

TYPE_DICT = {1 : vtk.VTK_LINE, 2 : vtk.VTK_TRIANGLE, 4 : vtk.VTK_TETRA,
             15 : vtk.VTK_PIXEL}

WRITER = {vtk.VTK_UNSTRUCTURED_GRID:vtk.vtkXMLUnstructuredGridWriter,
          vtk.VTK_POLY_DATA:vtk.vtkXMLPolyDataWriter,}

class GmshMesh(object):
    """This is a class for storing nodes and elements.

    Members:
    nodes -- A dict of the form { nodeID: [ xcoord, ycoord, zcoord] }
    elements -- A dict of the form { elemID: (type, [tags], [nodeIDs]) }

    Methods:
    read(file) -- Parse a Gmsh version 1.0 or 2.0 mesh file
    write(file) -- Output a Gmsh version 2.0 mesh file
    """

    def __init__(self):
        self.nodes = {}
        self.elements = {}

    def read(self, filename):
        """Read a Gmsh .msh file.

        Reads Gmsh format 1.0 and 2.0 mesh files, storing the nodes and
        elements in the appropriate dicts.
        """

        mshfile = open(filename, 'r')

        mode_dict = {'$NOD' : 1, '$Nodes' : 1,
                     '$ELM' : 2,
                     '$Elements' : 3}

        readmode = 0
        for line in mshfile:
            line = line.strip()
            if line.startswith('$'):
                readmode = mode_dict.get(line, 0)
            elif readmode:
                columns = line.split()
                if readmode == 1 and len(columns) == 4:
                    # Version 1.0 or 2.0 Nodes
                    try:
                        self.nodes[int(columns[0])] = [float(c) for c in columns[1:]]
                    except ValueError:
                        readmode = 0
                elif readmode > 1 and len(columns) > 5:
                    # Version 1.0 or 2.0 Elements
                    try:
                        columns = [int(c) for c in columns]
                    except ValueError:
                        readmode = 0
                    else:
                        (ele_id, ele_type) = columns[0:2]
                        if readmode == 2:
                            # Version 1.0 Elements
                            tags = columns[2:4]
                            nodes = columns[5:]
                        else:
                            # Version 2.0 Elements
                            ntags = columns[2]
                            tags = columns[3:3+ntags]
                            nodes = columns[3+ntags:]
                        self.elements[ele_id] = (ele_type, tags, nodes)

        mshfile.close()

    def write(self, filename):
        """Dump the mesh out to a Gmsh 2.0 msh file."""

        mshfile = open(filename, 'w')

        print >>mshfile, '$MeshFormat\n2.0 0 8\n$EndMeshFormat'
        print >>mshfile, '$Nodes\n%d'%len(self.nodes)
        for node_id, coord in self.nodes.items():
            print >>mshfile, node_id, ' '.join([str(c) for c in  coord])
        print >>mshfile, '$EndNodes'
        print >>mshfile, '$Elements\n%d'%len(self.elements)
        for ele_id, elem in self.elements.items():
            (ele_type, tags, nodes) = elem
            print >>mshfile, ele_id, ele_type, len(tags)
            print >>mshfile, ' '.join([str(c) for c in tags])
            print >>mshfile, ' '.join([str(c) for c in nodes])
        print >>mshfile, '$EndElements'

class PolyData(object):
    """ Class storing a living vtkPolyData construction"""

    def __init__(self, filename):
        """ Initialize the PolyData instance"""

        self.filename = filename
        self.cell_ids = {}
        self.poly_data = vtk.vtkPolyData()
        self.pnts = vtk.vtkPoints()
        self.pnts.Allocate(0)

        self.arrays = {}
        for name in ('Time',):
            array = vtk.vtkDoubleArray()
            array.SetName(name)
            self.poly_data.GetPointData().AddArray(array)

    def append_data(self, bucket):
        """ Add the data from the current time level"""

        for particle in bucket.particles:
            ids = self.cell_ids.setdefault(particle, vtk.vtkIdList())

            part_id = self.pnts.InsertNextPoint(particle.pos)

            ids.InsertNextId(part_id)

            self.poly_data.GetPointData().GetScalars('Time').InsertNextValue(bucket.time)

    def write(self):
        """ Write the staged vtkPolyData to a file."""

        self.poly_data.SetPoints(self.pnts)

        self.poly_data.Allocate(len(self.cell_ids))
        for cell_id in self.cell_ids.values():
            self.poly_data.InsertNextCell(vtk.VTK_LINE, cell_id)

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(self.filename)
        writer.SetInput(self.poly_data)
        writer.Write()

class BoundaryData(object):
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

    if any([celltypes.IsType(ct) for ct in TYPES_3D]):
        dim = 3
    elif any([celltypes.IsType(ct) for ct in TYPES_2D]):
        dim = 2
    elif any([celltypes.IsType(ct) for ct in TYPES_1D]):
        dim = 1
    else:
        dim = 0

    print dim

    ncells = ugrid.GetNumberOfCells()
    ncda = ugrid.GetCellData().GetNumberOfArrays()

    for i in range(ncda):
        out_grid.GetCellData().GetArray(i).SetName(ugrid.GetCellData().GetArray(i).GetName())

    cell_data = ugrid.GetCellData()
    for i in range(ncells):
        cell = ugrid.GetCell(i)
        if dim > cell.GetCellDimension():
            out_grid.InsertNextCell(cell.GetCellType(),
                                    cell.GetPointIds())
            for j in range(ncda):
                out_data = out_grid.GetCellData().GetArray(j)
                out_data.InsertNextTuple(cell_data.GetArray(j).GetTuple(i))

    return out_grid

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
    pos_x = []
    pos_y = []
    pos_z = []
    vel_u = []
    vel_v = []
    vel_w = []

    for line in infile.readlines():
        data = [float(M) for M in line.split()]

        time.append(data[0])
        num = len(data[1:])/6
        pos_x.append(data[1::3][:num])
        vel_u.append(data[1::3][num:])
        pos_y.append(data[2::3][:num])
        vel_v.append(data[2::3][num:])
        pos_z.append(data[3::3][:num])
        vel_w.append(data[3::3][num:])

    pos_x = numpy.array(pos_x)
    pos_y = numpy.array(pos_y)
    pos_z = numpy.array(pos_z)
    vel_u = numpy.array(vel_u)
    vel_v = numpy.array(vel_v)
    vel_w = numpy.array(vel_w)

    infile.close()

    return time, pos_x, pos_y, pos_z, vel_u, vel_v, vel_w


def ascii_to_polydata_time_series(filename, basename):

    """Convert ascii file to a series of vtkPolyData (.vtp) files.

    Each file contains one time level of the data, and are numbered sequentially.
    Within each file, each dataset is written to seperate pixel.

    Args:
        filename (str): Filename/path of the ascii file containing the data.
        basename (str): String used in the construction of the file series.
        The formula is of the form basename_0.vtp, basename_1.vtp,..."""

    ascii_data = get_ascii_data(filename)
    time = ascii_data[0]

    for i, full_data in enumerate(zip(ascii_data[1], ascii_data[2], ascii_data[3],
                                      ascii_data[4], ascii_data[5], ascii_data[6])):
        poly_data = vtk.vtkPolyData()
        pnts = vtk.vtkPoints()
        pnts.Allocate(0)
        poly_data.SetPoints(pnts)
        poly_data.Allocate(ascii_data[1].shape[0])

        outtime = vtk.vtkDoubleArray()
        outtime.Allocate(ascii_data[1].shape[0])
        outtime.SetName('Time')

        velocity = vtk.vtkDoubleArray()
        velocity.SetNumberOfComponents(3)
        velocity.Allocate(ascii_data[1].shape[0])
        velocity.SetName('Particle Velocity')

        for data in numpy.array(full_data).T:
            pixel = vtk.vtkPixel()
            pixel.GetPointIds().InsertId(0,
                                         poly_data.GetPoints().InsertNextPoint(data[0],
                                                                               data[1],
                                                                               data[2]))
            outtime.InsertNextValue(time[i])
            velocity.InsertNextTuple3(data[3], data[4], data[5])
            poly_data.InsertNextCell(pixel.GetCellType(), pixel.GetPointIds())

        poly_data.GetPointData().AddArray(outtime)
        poly_data.GetPointData().AddArray(velocity)
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName("%s_%d.vtp"%(basename, i))
        writer.SetInput(poly_data)
        writer.Write()

def write_level_to_polydata(bucket, level, basename, **kwargs):

    """Output a time level of a particle bucket to a vtkPolyData (.vtp) files.

    Each file contains one time level of the data, and are numbered sequentially.
    Within each file, each particle is written to seperate pixel.

    Args:
         bucket   (ParticleBucket):
        level    (int):
        basename (str): String used in the construction of the file series.
        The formula is of the form basename_0.vtp, basename_1.vtp,..."""

    del kwargs

    poly_data = vtk.vtkPolyData()
    pnts = vtk.vtkPoints()
    pnts.Allocate(0)
    poly_data.SetPoints(pnts)
    poly_data.Allocate(bucket.pos.shape[0])

    outtime = vtk.vtkDoubleArray()
    outtime.SetName('Time')
    outtime.Allocate(bucket.pos.shape[0])

    velocity = vtk.vtkDoubleArray()
    velocity.SetNumberOfComponents(3)
    velocity.Allocate(bucket.pos.shape[0])
    velocity.SetName('Particle Velocity')

    for positions, vel in zip(bucket.pos, bucket.vel):
        pixel = vtk.vtkPixel()
        pixel.GetPointIds().InsertId(0,
                                     poly_data.GetPoints().InsertNextPoint(positions[0],
                                                                           positions[1],
                                                                           positions[2]))
        outtime.InsertNextValue(bucket.time)
        velocity.InsertNextTuple3(vel[0], vel[1], vel[2])
        poly_data.InsertNextCell(pixel.GetCellType(), pixel.GetPointIds())

    poly_data.GetPointData().AddArray(outtime)
    poly_data.GetPointData().AddArray(velocity)
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName("%s_%d.vtp"%(basename, level))
    writer.SetInput(poly_data)
    writer.Write()


def write_level_to_ugrid(bucket, level, basename, model, **kwargs):
    """ Output a time level of a bucket to a vtkXMLUnstructuredGrid (.vtu) file.

    Particle volumes are averaged over the cells in the model using a control volume
    approach."""

    del kwargs

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.DeepCopy(model)

#    mass = numpy.zeros(ugrid.GetNumberOfPoints())
#    for _ in range(ugrid.GetNumberOfCells()):
#        cell = ugrid.GetCell(_)
#        dummy = cell.GetPoints().GetPoint
#        dummy_id = cell.GetPointIds().GetId
#
#        measure = abs(cell.TriangleArea(dummy(0), dummy(1), dummy(2)))
#
#        for k in range(cell.GetNumberOfPoints()):
#            mass[dummy_id(k)] += measure/cell.GetNumberOfPoints()

    locator = vtk.vtkPointLocator()
    locator.SetDataSet(ugrid)
    locator.BuildLocator()

    volume = numpy.zeros(ugrid.GetNumberOfPoints())
    temperature = numpy.zeros(ugrid.GetNumberOfPoints())
    velocity = numpy.zeros((ugrid.GetNumberOfPoints(), 3))

    LENGTH=0.1

    for dummy_particle in bucket.particles:
        point_list = vtk.vtkIdList()
        locator.FindPointsWithinRadius(LENGTH, dummy_particle.pos, point_list)

        for _ in range(point_list.GetNumberOfIds()):
            point_index = point_list.GetId(_)

            rad2 = numpy.sum((numpy.array(ugrid.GetPoints().GetPoint(point_index))
                              -dummy_particle.pos)**2)
            rad2 /= LENGTH

            volume[point_index] += (1.0/6.0*numpy.pi*dummy_particle.parameters.diameter**3
                                    *numpy.exp(-rad2**2))
            velocity[point_index, :] += (dummy_particle.vel*1.0/6.0*numpy.pi
                                         *dummy_particle.parameters.diameter**3
                                         *numpy.exp(-rad2**2))
        
    volume /= 0.5*LENGTH**2*(1.0-numpy.exp(-1.0**2))
    velocity /= 0.5*LENGTH**2*(1.0-numpy.exp(-1.0**2))

    for dummy_particle in bucket.particles:
        point_list = vtk.vtkIdList()
        locator.FindPointsWithinRadius(LENGTH, dummy_particle.pos, point_list)

        rad2 = numpy.sum((numpy.array(ugrid.GetPoints().GetPoint(point_index))
                              -dummy_particle.pos)**2)
        rad2 /= LENGTH

        for _ in range(point_list.GetNumberOfIds()):
            point_index = point_list.GetId(_)
            c = numpy.sum((dummy_particle.vel-velocity[point_index, :])**2)

        
            temperature [point_index] += (c*1.0/6.0*numpy.pi
                                         *dummy_particle.parameters.diameter**3
                                         *numpy.exp(-rad2**2))

    data = [vtk.vtkDoubleArray()]
    data[0].SetName('SolidVolumeFraction')
    data.append(vtk.vtkDoubleArray())
    data[1].SetName('SolidVolumeVelocity')
    data[1].SetNumberOfComponents(3)
    data.append(vtk.vtkDoubleArray())
    data[2].SetName('GranularTemperature')

    for _ in range(ugrid.GetNumberOfPoints()):
        data[0].InsertNextValue(volume[_])
        if volume[_] > 1.0e-20:
            data[1].InsertNextTuple3(*(velocity[_]/volume[_]))
        else:
            data[1].InsertNextTuple3(*(0*velocity[_]))
        data[2].InsertNextValue(temperature[_])
        ugrid.GetPointData().GetScalars('Time').SetValue(_, bucket.time)

    for _ in data:
        ugrid.GetPointData().AddArray(_)

    write_to_file(ugrid, "%s_out_%d.vtu"%(basename, level))

def ascii_to_polydata(filename, outfile):
    """Convert ascii file to a single vtkXMLPolyData (.vtp) files.

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

    for positions in zip(pos_x.T, pos_y.T, pos_z.T):
        line = vtk.vtkLine()
        for k, data in enumerate(zip(time, positions[0], positions[1], positions[2])):
            outtime.InsertNextValue(data[0])
            line.GetPointIds().InsertId(k,
                                        poly_data.GetPoints().InsertNextPoint(data[1],
                                                                              data[2],
                                                                              data[3]))
            poly_data.InsertNextCell(line.GetCellType(), line.GetPointIds())

    poly_data.GetPointData().AddArray(outtime)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(outfile)
    writer.SetInput(poly_data)

    writer.Write()

def collision_list_to_polydata(col_list, outfile,
                               model=Collision.mclaury_mass_coeff, **kwargs):
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
    normal = vtk.vtkDoubleArray()
    normal.SetNumberOfComponents(3)
    normal.SetName('Normal')

    for col in col_list:
        pixel = vtk.vtkPixel()
        pixel.GetPointIds().InsertId(0,
                                     poly_data.GetPoints().InsertNextPoint(col.pos[0],
                                                                           col.pos[1],
                                                                           col.pos[2]))
        time.InsertNextValue(col.time)
        wear.InsertNextValue(model(col, **kwargs))
        normal.InsertNextTuple3(*col.normal)
        poly_data.InsertNextCell(pixel.GetCellType(), pixel.GetPointIds())

    poly_data.GetPointData().AddArray(time)
    poly_data.GetPointData().AddArray(wear)
    poly_data.GetPointData().AddArray(normal)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(outfile)
    writer.SetInput(poly_data)

    writer.Write()

def get_linear_cell(cell):
    """ Get equivalent linear cell to vtkCell cell"""
    if cell.GetCellType() in (vtk.VTK_POLY_LINE,):
        linear_cell = vtk.vtkLine()
        linear_cell.GetPoints().SetPoint(0, cell.GetPoints().GetPoint(0))
        linear_cell.GetPoints().SetPoint(1, cell.GetPoints().GetPoint(1))
    elif cell.GetCellType() in (vtk.VTK_QUADRATIC_TRIANGLE,):
        linear_cell = vtk.vtkTriangle()
        linear_cell.GetPoints().SetPoint(0, cell.GetPoints().GetPoint(0))
        linear_cell.GetPoints().SetPoint(1, cell.GetPoints().GetPoint(1))
        linear_cell.GetPoints().SetPoint(2, cell.GetPoints().GetPoint(2))
    elif cell.GetCellType() in (vtk.VTK_QUADRATIC_TETRA,):
        linear_cell = vtk.vtkTetra()
        linear_cell.GetPoints().SetPoint(0, cell.GetPoints().GetPoint(0))
        linear_cell.GetPoints().SetPoint(1, cell.GetPoints().GetPoint(1))
        linear_cell.GetPoints().SetPoint(2, cell.GetPoints().GetPoint(2))
        linear_cell.GetPoints().SetPoint(3, cell.GetPoints().GetPoint(3))
    else:
        linear_cell = cell

    return linear_cell

def test_in_cell(cell, position):
    """ Check if point is in vtk cell"""

    linear_cell = get_linear_cell(cell)
    dim = linear_cell.GetNumberOfPoints()-1
    ppos = numpy.zeros(linear_cell.GetNumberOfPoints())
    dummy_func = cell.GetPoints().GetPoint
    args = [dummy_func(i)[:dim] for i in range(1, dim+1)]
    args.append(dummy_func(0)[:dim])
    args.append(ppos)
    linear_cell.BarycentricCoords(position[:dim], *args)

    return cell.GetParametricDistance(ppos[:3]) == 0

def write_to_file(vtk_data, outfile):
    """ Wrapper around the various VTK writer routines"""

    writer = WRITER[vtk_data.GetDataObjectType()]()
    writer.SetFileName(outfile)
    writer.SetInput(vtk_data)
    writer.Write()

def make_unstructured_grid(mesh, velocity, pressure, time, outfile=None):
    """Given a mesh (in Gmsh format), velocity and pressure fields, and a
    time level, store the data in a vtkUnstructuredGridFormat."""

    pnts = vtk.vtkPoints()
    pnts.Allocate(len(mesh.nodes))

    node2id = {}
    for k, point in mesh.nodes.items():
        node2id[k] = pnts.InsertNextPoint(point)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(pnts)

    for element in mesh.elements.values():
        id_list = vtk.vtkIdList()
        for node in element[2]:
            id_list.InsertNextId(node2id[node])

        ugrid.InsertNextCell(TYPE_DICT[element[0]], id_list)

    data = []
    data.append(vtk.vtkDoubleArray())
    data[0].SetNumberOfComponents(3)
    data[0].Allocate(3*pnts.GetNumberOfPoints())
    data[0].SetName('Velocity')

    data.append(vtk.vtkDoubleArray())
    data[1].Allocate(pnts.GetNumberOfPoints())
    data[1].SetName('Pressure')

    data.append(vtk.vtkDoubleArray())
    data[2].Allocate(pnts.GetNumberOfPoints())
    data[2].SetName('Time')

    for k in range(len(mesh.nodes)):
        if hasattr(velocity, '__call__'):
            data[0].InsertNextTuple3(*velocity(ugrid.GetPoints().GetPoint(k)))
        else:
            data[0].InsertNextTuple3(*velocity[k, :])
        if hasattr(pressure, '__call__'):
            data[1].InsertNextValue(pressure(ugrid.GetPoints().GetPoint(k)))
        else:
            data[1].InsertNextValue(pressure[k])
        data[2].InsertNextValue(time)


    for _ in data:
        ugrid.GetPointData().AddArray(_)

    if outfile:
        write_to_file(ugrid, outfile)

    return ugrid

def interpolate_collision_data(col_list, ugrid, method='nearest'):
    """ Interpolate the wear data from col_list onto the surface described
    in the vtkUnstructuredGrid ugrid """

    out_pts = numpy.array([ugrid.GetPoint(i) for i in range(ugrid.GetNumberOfPoints())])

    data_pts = numpy.array([c.pos for c in col_list])
    wear_pts = numpy.array([c.get_wear() for c in col_list])

    print data_pts, numpy.isfinite(data_pts).all()
    print wear_pts, numpy.isfinite(wear_pts).all()

    wear_on_grid = griddata(data_pts, wear_pts, out_pts, method=method)
#    rbfi = Rbf(data_pts[:,0], data_pts[:, 1], data_pts[:, 2], wear_pts,
#    function='gaussian', epsilon=0.1,smooth=3.0)
#    wear_on_grid = rbfi(out_pts[:,0], out_pts[:, 1], out_pts[:, 2])

    print wear_on_grid.max()

    wear_vtk = vtk.vtkDoubleArray()
    wear_vtk.SetName('Wear')

    for wear in wear_on_grid:
        wear_vtk.InsertNextValue(wear)

    ugrid.GetPointData().AddArray(wear_vtk)
