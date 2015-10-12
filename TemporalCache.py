""" Module controlling access to the vtu/pvtu files driving the particle model"""

import vtk
import glob
import numpy

class TemporalCache(object):
    """ The base object containing the vtu files. 

    This scans the files for their timelevel information and provides a pair of files bracketing the desired timelevel when called
    """
    def __init__(self,base_name,t_min=0.,t_max=numpy.infty):
        """ 
        Initialise the cache from a base file name and optional limits on the time levels desired.
        """
        files=glob.glob(base_name+'*.vtu')

        self.data = []
        self.reset()

        print files

        for filename in files:
            rdr = vtk.vtkXMLUnstructuredGridReader()
            rdr.SetFileName(filename)
            for k in range(rdr.GetNumberOfPointArrays()):
                rdr.SetPointArrayStatus(rdr.GetPointArrayName(k), 0)
            for k in range(rdr.GetNumberOfCellArrays()):
                rdr.SetCellArrayStatus(rdr.GetCellArrayName(k), 0)
            rdr.SetPointArrayStatus('Time', 1)
            rdr.Update()
            ug = rdr.GetOutput()
            t = ug.GetPointData().GetScalars('Time').GetValue(0)
            
            self.data.append([t, filename,None,None])

        self.data.sort(cmp=lambda x, y:cmp(x[0], y[0]))

        self.range(t_min, t_max)
            
    def reset(self):
        """ Reset the bounds on the loaded cache data"""
        for k,d in enumerate(self.data):
            if d[2]: self.close(k)
        self.lower = 0
        self.upper = 0

    def range(self,a,b):
        """ Specify a range of data to keep open."""
        if len(self.data) == 0: return
        if self.data[self.lower][0] > a:
            self.reset()
        while self.lower < len(self.data)-2 and self.data[self.lower+1][0] <= a:
            self.close(self.lower)
            self.lower += 1
        if self.upper <= self.lower: 
            self.upper = self.lower
            self.open(self. lower)
        while self.upper <= len(self.data)-2 and self.data[self.upper][0] <= b:
            self.upper += 1
            self.open(self.upper)

        return self.data[self.lower:self.upper+1]

    def open(self, k):
        """ Open a file for reading."""
        rdr=vtk.vtkXMLUnstructuredGridReader()

        print 'loading %s'%self.data[k][1]
        rdr.SetFileName(self.data[k][1])
        rdr.Update()

        self.data[k][2] = rdr.GetOutput()
        cl=vtk.vtkCellLocator()
        cl.SetDataSet(self.data[k][2])
        cl.BuildLocator()
        self.data[k][3] = cl

    def close(self, k):
        """Close an open file (implictly through the garbage collector."""
        del self.data[k][3]
        del self.data[k][2]
        self.data[k].append(None)
        self.data[k].append(None)

    def __call__(self, t):
        """ Get the data bracketing time level t."""
        lower = self.lower
        upper = self.upper

        assert self.data[lower][0] <= t and self.data[upper][0] >= t
        
        while lower < len(self.data)-2 and self.data[lower+1][0] <= t:
            lower+=1

        t0 = self.data[lower][0]
        t1 = self.data[lower+1][0]
        if t1 == t0: t1 = numpy.infty

        return self.data[lower:lower+2], (t-t0)/(t1-t0)
