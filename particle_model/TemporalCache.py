""" Module controlling access to the vtu/pvtu files driving the particle model"""

import glob
import vtk
import numpy
from particle_model import Parallel
from particle_model import Debug
from Debug import profile
from particle_model import vtk_extras
try:
    from lxml import etree as ET
    def element_tree(**kwargs):
        """ Wrapper for lxml ElementTree."""
        return ET.ElementTree(parser=ET.XMLParser(recover=True), **kwargs)
except ImportError:
    from xml.etree import ElementTree as ET
    def element_tree(**kwargs):
        """ Wrapper for xml ElementTree."""
        return ET.ElementTree(**kwargs)

PICKERS = [ vtk_extras.Picker(),
            vtk_extras.Picker(),
            vtk_extras.Picker()]

def read_pvd(filename):
    """Read timestep and filename data from a .pvd file."""
    times = []
    names = []
    etree = element_tree(file=filename).getroot()

    for data in etree[0].findall('DataSet'):
        times.append(float(data.get('timestep')))
        names.append((data.get('file')))

    return zip(times, names)

def get_piece_filename_from_vtk(filename, piece=Parallel.get_rank()):

    """Get the filename of individual VTK file piece."""

    etree = element_tree(file=filename).getroot()
    return etree[0].findall('Piece')[piece].get('Source')

class DataCache(object):
    """ Store data in a cyclical cache. """

    def __init__(self, max_size=2):
        self._keys = [None for _ in range(max_size)]
        self._values = [{} for _ in range(max_size)]

    @profile
    def get(self, infile, name):
        """Get cache value."""
        try:
            _ = self._values[self._keys.index(infile)]
        except ValueError:
            self._keys[1:] = self._keys[:-1]
            self._values[1:] = self._values[:-1]
            self._keys[0] = infile
            self._values[0] = {}
            _ = self._values[0]

        data = _.get(name)
        if not data:
            if infile.IsA('vtkUnstructuredGrid'):
                if infile.GetPointData().HasArray(name):
                    data = infile.GetPointData().GetArray(name)
                else:
                    data = None
            else:
                data = None
                for _ in range(infile.GetNumberOfBlocks()):
                    if infile.GetBlock(_).GetPointData().HasArray(name):
                        data = infile.GetBlock(_).GetPointData().GetArray(name)
                        break
            self.set(infile, name, data)

        return data

    def set(self, infile, name, data):
        """Set cache value."""
        try:
            _ = self._values[self._keys.index(infile)]
        except ValueError:
            self._keys[1:] = self._keys[:-1]
            self._values[1:] = self._values[:-1]
            self._keys[0] = infile
            self._values[0] = {}
            _ = self._values[0]

        _[name] = data

class TemporalCache(object):
    """ The base object containing the vtu files.

    This scans the files for their timelevel information and provides a pair
    of files bracketing the desired timelevel when called
    """
    def __init__(self, base_name, t_min=0., t_max=numpy.infty, online=False,
                 parallel_files=False, timescale_factor=1.0, **kwargs):
        """
        Initialise the cache from a base file name and optional limits on the time levels desired.
        """

        self.data = []
        self.set_field_names(**kwargs)
        self.reset()

        if base_name.rsplit(".", 1)[-1] == "pvd":
            for time, filename in read_pvd(base_name):
                self.data.append([timescale_factor*time, filename, None, None])
        else:
            if (Parallel.is_parallel() and online) or parallel_files:
                files = glob.glob(base_name+'_[0-9]*.pvtu')
            else:
                files = glob.glob(base_name+'_[0-9]*.vtu')

            for filename in files:
                if (Parallel.is_parallel() and online) or parallel_files:
                    pfilename = get_piece_filename_from_vtk(filename)
                else:
                    pfilename = filename
                time = self.get_time_from_vtk(pfilename)
                self.data.append([timescale_factor*time, pfilename, None, None])

        self.data.sort(cmp=lambda x, y: cmp(x[0], y[0]))
        self.range(t_min, t_max)

        self.cache = DataCache()

    def set_field_names(self, velocity_name="Velocity", pressure_name="Pressure", time_name="Time"):
        """Set the names used to look up pressure and velocity fields."""
        self.field_names = {}
        self.field_names["Velocity"] = velocity_name or ""
        self.field_names["Pressure"] = pressure_name or ""
        self.field_names["Time"] = time_name or ""

    def get(self, infile, name):
        """Find array, possibly from cache."""
        return self.cache.get(infile, name)

    def reset(self):
        """ Reset the bounds on the loaded cache data"""
        for k, dat in enumerate(self.data):
            if dat[2]:
                self.close(k)
        self.lower = 0
        self.upper = 0

    def range(self, t_min, t_max):
        """ Specify a range of data to keep open."""
        if not self.data:
            return
        if self.data[self.lower][0] > t_min:
            self.reset()
        while (self.lower < len(self.data)-2
               and self.data[self.lower+1][0] <= t_min):
            self.close(self.lower)
            self.lower += 1
        if self.upper <= self.lower:
            self.upper = self.lower
            self.open(self. lower)
        while (self.upper <= len(self.data)-2
               and self.data[self.upper][0] <= t_max):
            self.upper += 1
            self.open(self.upper)

        return self.data[self.lower:self.upper+1]

    def open(self, k):
        """ Open a file for reading."""
        rdr = vtk.vtkXMLGenericDataObjectReader()

        Debug.logger.info('loading %s'%self.data[k][1])
        rdr.SetFileName(self.data[k][1])
        rdr.Update()

        self.data[k][2] = rdr.GetOutput()
        cloc = vtk.vtkCellLocator()
        cloc.SetDataSet(self.data[k][2])
        cloc.BuildLocator()
        self.data[k][3] = cloc

    def close(self, k):
        """Close an open file (implictly through the garbage collector."""
        del self.data[k][3]
        del self.data[k][2]
        self.data[k].append(None)
        self.data[k].append(None)

    def get_time_from_vtk(self, filename):
        """ Get the time from a vtk XML formatted file."""

        parallel_files = ('pvtu', 'pvtp', 'pvtm', 'vtm')
        parallel = filename.split('.')[-1] in parallel_files

        etree = element_tree(file=filename).getroot()
        assert etree.tag == 'VTKFile'
        if parallel:
            parallel_name = etree[0].findall('Piece')[Parallel.get_rank()].get('Source')
            return self.get_time_from_vtk(parallel_name)
        else:
            for piece in etree[0]:
                for data in piece[0]:
                    if data.get('Name') != self.field_names["Time"]:
                        continue
                    return float(data.get('RangeMin'))

    def __call__(self, time):
        """ Get the data bracketing time level."""
        lower = self.lower
        upper = self.upper

        assert self.data[lower][0] <= time and self.data[upper][0]+1.0e-8 >= time

        while lower < len(self.data)-2 and self.data[lower+1][0] <= time:
            lower += 1

        t_min = self.data[lower][0]
        t_max = self.data[lower+1][0]
        if t_max == t_min:
            t_max = numpy.infty

        return (self.data[lower:lower+2], (time-t_min)/(t_max-t_min),
                [[self.field_names["Velocity"], self.field_names["Pressure"]],
                 [self.field_names["Velocity"], self.field_names["Pressure"]]])

    def __iter__(self):
        return self.data.__iter__()

    def get_bounds(self, ptime):
        """ Get bounds of vtk object, in form (xmin, xmax, ymin, ymax, zmin, zmax)."""
        data = self(ptime)
        bounds = numpy.zeros(6)
        data[0][1][-2].ComputeBounds()
        data[0][1][-2].GetBounds(bounds)
        return bounds

    def get_velocity(self, pos, time):
        """ Get the velocity value from the cache at a given position and time."""

        data, alpha, names = self(time)

        vel_data0 = self.cache.get(data[0][2],names[0][0])
        vel_data1 = self.cache.get(data[1][2],names[1][0])

        loc0 = data[0][3]
        loc0.BuildLocatorIfNeeded()
        loc1 = data[1][3]
        loc1.BuildLocatorIfNeeded()

        vel0 = vtk_extras.EvaluateField(vel_data0, loc0, pos, names[0][0])
        vel1 = vtk_extras.EvaluateField(vel_data1, loc1, pos, names[1][0])

        return alpha*vel1+(1.0-alpha)*vel0

class FluidityCache(object):
    """Cache like object used when running particles online."""

    def __init__(self, block, time, dt, velocity_name='Velocity'):
        """ Initialise the cache.
             block  -- The VTK multiblock object
             time   -- The (current) simulation time
             dt     -- The model timestep """
        self.block = block
        self.time = time
        self.delta_t = dt
        self.velocity_name = velocity_name
        self.cloc = vtk.vtkCellLocator()
        self.cloc.SetDataSet(self.block.GetBlock(0))
        self.cloc.SetTolerance(0.0)
        self.cloc.BuildLocator()
        self.cache = DataCache()

    def get(self, infile, name):
        """Get array, possibly from cache."""
        return self.cache.get(infile, name)

    def update(self, block, time, delta_t):
        """Update latest time level of a fluidity style data cache."""
        self.block = block
        self.time = time
        self.delta_t = delta_t
        self.cloc.SetDataSet(self.block.GetBlock(0))
        self.cloc.SetTolerance(0.0)
        self.cloc.BuildLocator()

    def range(self, t_min, t_max):
        """ Specify a range of data to keep open. Not used here"""
        pass

    def __call__(self, ptime):
        self.cloc.Update()
        return ([[self.time, None, self.block, self.cloc],
                 [self.time+self.delta_t, None, self.block, self.cloc]],
                (ptime-self.time+self.delta_t)/self.delta_t,
                [['Old'+self.velocity_name, 'OldPressure'],
                 [self.velocity_name, 'Pressure']])

    def get_bounds(self, ptime):
        """ Get the bounds of vtk object,  in form (xmin, xmax, ymin, ymax, zmin, zmax)."""
        del ptime
        bounds = numpy.zeros(6)
        self.block.GetBlock(0).ComputeBounds()
        self.block.GetBlock(0).GetBounds(bounds)
        return bounds
