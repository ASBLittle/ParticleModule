import vtk
import glob
import numpy

class TemporalCache(object):
    
    def __init__(self,base_name,t_min=0.,t_max=numpy.infty):
        files=glob.glob(base_name+'*.vtu')

        self.data=[]
        self.reset()

        print files

        for file in files:
            rdr=vtk.vtkXMLUnstructuredGridReader()
            rdr.SetFileName(file)
            for k in range(rdr.GetNumberOfPointArrays()):
                rdr.SetPointArrayStatus(rdr.GetPointArrayName(k),0)
            for k in range(rdr.GetNumberOfCellArrays()):
                rdr.SetCellArrayStatus(rdr.GetCellArrayName(k),0)
            rdr.SetPointArrayStatus('Time',1)
            rdr.Update()
            ug=rdr.GetOutput()
            t=ug.GetPointData().GetScalars('Time').GetValue(0)
            
            self.data.append([t,file,None,None])

        self.data.sort(cmp=lambda x,y:cmp(x[0],y[0]))

        self.range(t_min,t_max)
            
    def reset(self):
        for k,d in enumerate(self.data):
            if d[2]: self.close(k)
        self.lower=0
        self.upper=0

    def range(self,a,b):
        if len(self.data)==0: return
        if self.data[self.lower][0]>a:
            self.reset()
        while self.lower<len(self.data)-2 and self.data[self.lower+1][0]<=a:
            self.close(self.lower)
            self.lower+=1
        if self.upper<=self.lower: 
            self.upper=self.lower
            self.open(self.lower)
        while self.upper<=len(self.data)-2 and self.data[self.upper][0]<=b:
            self.upper+=1
            self.open(self.upper)

        return self.data[self.lower:self.upper+1]

    def open(self,k):
        
        rdr=vtk.vtkXMLUnstructuredGridReader()

        print 'loading %s'%self.data[k][1]
        rdr.SetFileName(self.data[k][1])
        rdr.Update()

        self.data[k][2]=rdr.GetOutput()
        cl=vtk.vtkCellLocator()
        cl.SetDataSet(self.data[k][2])
        cl.BuildLocator()
        self.data[k][3]=cl

    def close(self,k):
        del self.data[k][3]
        del self.data[k][2]
        self.data[k].append(None)
        self.data[k].append(None)

    def __call__(self,t):
        lower=self.lower
        upper=self.upper
        
        while lower<len(self.data)-2 and self.data[lower+1][0]<=t:
            lower+=1

        t0=self.data[lower][0]
        t1=self.data[lower+1][0]
        if t1==t0: t1=numpy.infty

        return self.data[lower:lower+2], (t-t0)/(t1-t0)
