""" Module controlling access to the vtu/pvtu files driving the particle model"""

import glob
import numpy
import vtk
from particle_model import Parallel
from particle_model import Debug
from particle_model import IO
from particle_model.Debug import profile
from particle_model import vtk_extras
from scipy.io import netcdf
import scipy
import utm
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

PICKERS = [vtk_extras.Picker(),
           vtk_extras.Picker(),
           vtk_extras.Picker()]

def map_bathy(bathy_file, var_name='z'):
    nc = netcdf.NetCDFFile(bathy_file,'r', mmap=False)
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    values = nc.variables[var_name][:,:]
    bathymetry_interpolator = scipy.interpolate.RegularGridInterpolator((lat, lon), values)
    nc.close()
    return bathymetry_interpolator

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
            if (infile.IsA('vtkUnstructuredGrid') or
                infile.IsA('vtkStructuredGrid') or
                infile.IsA('vtkRectilinearGrid')):
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
                files = glob.glob(base_name+'_[0-9]*.p%s'%kwargs.get('fileext',
                                                                     'vtu'))
            else:
                files = glob.glob(base_name+'_[0-9]*.%s'%kwargs.get('fileext',
                                                                     'vtu'))

            for filename in files:
                if (Parallel.is_parallel() and online) or parallel_files:
                    pfilename = get_piece_filename_from_vtk(filename)
                else:
                    pfilename = filename
                time = self.get_time_from_vtk(pfilename)
                self.data.append([timescale_factor*time, pfilename, None, None])

        self.data.sort(key=lambda x: x[0])
        self.range(t_min, t_max)

        self.cache = DataCache()

    def set_field_names(self, velocity_name="Velocity",
                        pressure_name="Pressure", time_name="Time",
                        **kwargs):
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
            raise ValueError
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

        Debug.logger.info('loading %s', self.data[k][1])
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

        parallel_files = ('pvtu', 'pvtp', 'pvtm', 'vtm', 'pvts', 'pvtr')
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

        vel_data0 = self.cache.get(data[0][2], names[0][0])
        vel_data1 = self.cache.get(data[1][2], names[1][0])

        loc0 = data[0][3]
        loc0.BuildLocatorIfNeeded()
        loc1 = data[1][3]
        loc1.BuildLocatorIfNeeded()

        vel0 = vtk_extras.EvaluateField(vel_data0, loc0, pos, names[0][0])
        vel1 = vtk_extras.EvaluateField(vel_data1, loc1, pos, names[1][0])

        return alpha*vel1+(1.0-alpha)*vel0


    def get_mod_velocity(self, pos, z, sigma, time, surface_depth, correction_factor=1, wind_VD=[0,0], prandle=[-0.61, 1.16, 0.62]):
        '''
        Modify base velocity to include vertical profile from depth averaged
        velocity value - using prandle profile
        (default is a modified profile as used in PC-GESTS)
        calculate wind driven flow in upper water
        z_c is depth below which wind does not influnce flow
        c_d is fractional wind speed tranfer factor default 0.035
        wind_VD - velocity and bearing as one
        tidegauge elev is used to determine free surface height - single value
        used to reduce function calls

        '''
        data, alpha, names = self(time)
        vel_data0 = self.cache.get(data[0][2], names[0][0])
        vel_data1 = self.cache.get(data[1][2], names[1][0])
        loc0 = data[0][3]
        loc0.BuildLocatorIfNeeded()
        loc1 = data[1][3]
        loc1.BuildLocatorIfNeeded()
        vel0 = vtk_extras.EvaluateField(vel_data0, loc0, pos, names[0][0])
        vel1 = vtk_extras.EvaluateField(vel_data1, loc1, pos, names[1][0])
        da_vel = alpha*vel1+(1.0-alpha)*vel0

        zeta = 1-sigma
        # zeta = 1 - sigdep(pos, dep, tidalelev)  # 1- this value if required
        # vertical velocity profile - in this usage zeta = 1 - sigma
        z_vel = da_vel * (prandle[0]*zeta**2 + prandle[1]*zeta + prandle[2])  # only works for arrays
        # correction of -2* z_vel_corr in gareloch above 55*, below 10m deep
        z_vel[2] = da_vel[2]  # reset vertical velocity incorrectly adjusted above
        # wind data
        c_d = 0.035  # SURFACE WIND CURRENT FACTOR
        z_0 = 0.10  + surface_depth  # SURFACE ROUGHNESS LENGTH (M) plus surface elevation depth (true z c.f. model 0)
        alpha = 2.00  # WAVE INFLUENCE FACTOR
        # FOR FUTURE - option to update wave data if desired
        # ASSUME A 3S WAVE IF NO WAVE DATA ARE PROVIDED
        t_z=3.  # wave period (s)
        r_l=9.81*t_z**2/(2.*numpy.pi)  # roughness length?
        z_c = (alpha * r_l) + surface_depth  # (true z c.f. model 0)
        # Calculate wind influence
        if wind_VD[0] == 0 or z > z_c:  # no wind or below influence depth
            U_wnd = 0
        else:
            U_s = c_d * wind_VD[0]  # velocity factor * wind speed (magnitude)
            if z < z_0:  # for surface levels only
                U_wnd = U_s
            else:
                U_wnd = U_s * (1-numpy.log(z/z_0)/numpy.log(z_c/z_0))
        z_vel[0] += U_wnd*numpy.sin(numpy.radians(wind_VD[1]-180))  # x
        z_vel[1] += U_wnd*numpy.cos(numpy.radians(wind_VD[1]-180))  # y
        return z_vel * correction_factor

class Elev_TempCache(object):
    """ The base object containing the vtu files.

    This scans the files for their timelevel information and provides a pair
    of files bracketing the desired timelevel when called
    """
    def __init__(self, base_name, bathymetry_interpolator, t_min=0., t_max=numpy.infty, online=False,
                  parallel_files=False, timescale_factor=1.0, minimum_depth = 10, utm_zone=30, utm_band='U', **kwargs):
        """
        Initialise the cache from a base file name and optional limits on the time levels desired.
        """
        self.minimum_depth = minimum_depth  # Fixed height above MSL to avoid negative depths
        self.bathymetry_interpolator = bathymetry_interpolator  # fixed bathymetry depth WRT MSL
        self.utm_zone=utm_zone 
        self.utm_band=utm_band
        self.data = []
        self.set_field_names(**kwargs)
        self.reset()

        if base_name.rsplit(".", 1)[-1] == "pvd":
            for time, filename in read_pvd(base_name):
                self.data.append([timescale_factor*time, filename, None, None])
        else:
            if (Parallel.is_parallel() and online) or parallel_files:
                files = glob.glob(base_name+'_[0-9]*.p%s'%kwargs.get('fileext',
                                                                     'vtu'))
            else:
                files = glob.glob(base_name+'_[0-9]*.%s'%kwargs.get('fileext',
                                                                     'vtu'))

            for filename in files:
                if (Parallel.is_parallel() and online) or parallel_files:
                    pfilename = get_piece_filename_from_vtk(filename)
                else:
                    pfilename = filename
                time = self.get_time_from_vtk(pfilename)
                self.data.append([timescale_factor*time, pfilename, None, None])

        self.data.sort(key=lambda x: x[0])
        self.range(t_min, t_max)

        self.cache = DataCache()

    def set_field_names(self, elevation_name="Elevation", time_name="Time",
                        **kwargs):
        """Set the names used to look up pressure and velocity fields."""
        self.field_names = {}
        self.field_names["Elevation"] = elevation_name or ""
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
            raise ValueError
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

        Debug.logger.info('loading %s', self.data[k][1])
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

        parallel_files = ('pvtu', 'pvtp', 'pvtm', 'vtm', 'pvts', 'pvtr')
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
                [[self.field_names["Elevation"]],
                 [self.field_names["Elevation"]]])

    def __iter__(self):
        return self.data.__iter__()

    def get_bounds(self, ptime):
        """ Get bounds of vtk object, in form (xmin, xmax, ymin, ymax, zmin, zmax)."""
        data = self(ptime)
        bounds = numpy.zeros(6)
        data[0][1][-2].ComputeBounds()
        data[0][1][-2].GetBounds(bounds)
        return bounds

    def get_elevation(self, pos, time):
        """ Get the free surface elevation from the cache at a given position and time.
        elevation is c.f. Mean Sea Level datum but returned value is c.f. model datum
        i.e. MSL + minimum_depth (default 10m) 
        ---- 0m
        ==== free surface
        <elevation>
        ---- MSL
        <bathy depth>
        === seabed
        """

        data, alpha, names = self(time)

        elev_data0 = self.cache.get(data[0][2], names[0][0])
        elev_data1 = self.cache.get(data[1][2], names[1][0])

        loc0 = data[0][3]
        loc0.BuildLocatorIfNeeded()
        loc1 = data[1][3]
        loc1.BuildLocatorIfNeeded()

        elev0 = vtk_extras.EvaluateField(elev_data0, loc0, pos, names[0][0])
        elev1 = vtk_extras.EvaluateField(elev_data1, loc1, pos, names[1][0])
        # interpolate between vtu timesteps
        elev = (alpha*elev1+(1.0-alpha)*elev0)[0]
        # return z and dep - wrt model zero - MSL + min_depth and MSL respectively
        return self.minimum_depth - elev, elev

    def get_depth(self, pos, time):
        """ Get the total depth from the cache at a given position and time.
        returns z - total depth and water column depth
        """
        # Bathymetry data at position
        bathy_dep = self.get_bathy(pos)
        # Elevation Data
        elev_z, elev_d = self.get_elevation(pos, time)
        # total depth in z, water depth in dep
        return bathy_dep + self.minimum_depth, bathy_dep + elev_d

    def get_sigma_depth(self, pos, z, time):
        """ Get the sigma depth value from the cache at a given position and time.
        Also returns water depth at the particle location as this is essentially
        free and could be useful later"""
        # Bathymetry data at position
        depth_z, water_depth = self.get_depth(pos, time)
        surface_depth = depth_z - water_depth
        part_depth = z - surface_depth # depth below free surface
        return max(min(1., part_depth / water_depth), 0.), surface_depth, water_depth  # limits return to between 0 and 1 even if true z is above or below min/max water depth
        # return z / dep  # unbounded to allow setting particle data

    def set_z_on_sigma(self, pos, sigma, time):
        depth_z, water_depth = self.get_depth(pos, time)
        surface_depth = depth_z - water_depth
        part_depth = sigma * water_depth  # c.f. MSL
        part_z = part_depth + surface_depth  # c.f. model 0
        return part_z, depth_z, surface_depth, water_depth

    def get_bathy(self, pos):
        try:
            lat, lon = utm.to_latlon(pos[0], pos[1], self.utm_zone, self.utm_band)
        except:
            return 0
        return -self.bathymetry_interpolator((lat, lon))  # negative depths in original


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

    def get_velocity(self, pos, time):
        """ Get the velocity value from the cache at a given position and time."""

        data, alpha, names = self(time)

        data[0][3].BuildLocatorIfNeeded()
        data[1][3].BuildLocatorIfNeeded()

        PICKERS[0].name = names[0][0]
        PICKERS[0].grid = IO.get_block(data[0][2], names[0][0])
        PICKERS[0].locator = data[0][3]
        PICKERS[0].pos = pos
        PICKERS[1].name = names[1][0]
        PICKERS[1].grid = IO.get_block(data[1][2], names[1][0])
        PICKERS[1].locator = data[1][3]
        PICKERS[1].pos = pos

        vel0 = PICKERS[0].nearest(pos)
        vel1 = PICKERS[1].nearest(pos)

        if vel0 is None or vel1 is None:
            raise ValueError
        #otherwise
        return alpha*vel1+(1.0-alpha)*vel0
