import vtk
import numpy
import pylab as p
import itertools


reader=vtk.vtkXMLUnstructuredGridReader()
bnd_reader=vtk.vtkXMLUnstructuredGridReader()

file1=vtk.vtkUnstructuredGrid()
file2=vtk.vtkUnstructuredGrid()
l1=vtk.vtkCellLocator()
l2=vtk.vtkCellLocator()

bnd_reader.SetFileName('bnd.vtu')
bfile=bnd_reader.GetOutput()
bnd_reader.GetOutput().Update()
bfile.Update()
gf=vtk.vtkGeometryFilter()
gf.SetInput(bfile)
gf.Update()
bndl=vtk.vtkCellLocator()
bndl.SetDataSet(gf.GetOutput())
bndl.BuildLocator()

class collisionException(Exception):
    pass

class collisionInfo(object):
    def __init__(self,x,cell,angle,t):
        self.x=x
        self.cell=cell
        self.angle=angle
        self.t=t

class particle(object):

    def __init__(self,p,v,t=0.0,dt=1.0):
        
        self.p=p
        self.v=v
        self.t=t
        self.dt=dt
        self.collisions=[]

    def update(self):
        # Simple RK4 integration

        collision=False

        try:
            kp1=self.v
            kv1=self.force(self.p,
                       self.v,
                       self.t)
        except:
            print self.p

        s,c1=self.collide(self.p+0.5*self.dt*kp1)
        kp2=self.v+0.5*self.dt*kv1
        kv2=self.force(self.p+s,
                       self.v+0.5*self.dt*kv1,
                       self.t+0.5*self.dt)

        s,c2=self.collide(self.p+0.5*self.dt*kp2)
        kp3=self.v+0.5*self.dt*kv2
        kv3=self.force(self.p+s,
                       self.v+0.5*self.dt*kv2,
                       self.t+0.5*self.dt)


        s,c3=self.collide(self.p+self.dt*kp3)
        kp4=self.v+self.dt*kv3
        kv4=self.force(self.p+s,
                       self.v+self.dt*kv3,
                       self.t+self.dt)


        s,c4,v4=self.collide(self.p+self.dt*(kp1+2.0*(kp2+kp3)+kp4)/6.0,
                             v=self.v+self.dt*(kv1+2.0*(kv2+kv3)+kv4)/6.0)
        self.p+=s
        self.v+=v4
        if type(c4)!=type(None):
            self.collisions.append(c4)


        self.t+=self.dt

    def force(self,p,v,t):
        
        ## Drag term:

        beta=0.44

        u,collision=self.picker(p,t)

#        if collision:
#            raise collisionException

        drag=beta*(u-v)

        return drag


    def picker(self,p,t):

        

        l1.BuildLocatorIfNeeded()
        l2.BuildLocatorIfNeeded()

        def fpick(file,locator):

            cp=[0.0,0.0,0.0]
            ci=vtk.mutable(0)
            si=vtk.mutable(0)
            sd=vtk.mutable(0.0)

            locator.FindClosestPoint(p,cp,ci,si,sd)

            c=file.GetCell(ci)

            collision=not c.PointInTriangle(p,c.GetPoints().GetPoint(0),
                              c.GetPoints().GetPoint(1),
                              c.GetPoints().GetPoint(2),1.0e-6)
                

            v1=[0.0,0.0]
            v2=[0.0,0.0]
            v3=[0.0,0.0]

            c.ProjectTo2D(c.GetPoints().GetPoint(0),
                          c.GetPoints().GetPoint(1),
                          c.GetPoints().GetPoint(2),
                          v1,v2,v3)

            x=[0.0,0.0,0.0]

            c.BarycentricCoords(p[:2],v1,v2,v3,x)
        

            pids=c.GetPointIds()
            
            data=file.GetPointData().GetVectors('Velocity')

            out=numpy.zeros(3)
            
            for k,a in enumerate(x):
                out+=a*numpy.array(data.GetTuple(pids.GetId(k)))

            return out,collision
            
        return fpick(file1,l1)


    def collide(self,p,v=None,pa=None,dt=None):
        if type(pa)==type(None):
            pa=self.p
            dt=self.dt
        e=0.9
        

        s=vtk.mutable(-1.0)
        x=[0.0,0.0,0.0]
        loc=[0.0,0.0,0.0]
        si=vtk.mutable(0)
        ci=vtk.mutable(0)

        f=bndl.IntersectWithLine(pa,p,
                               1.0e-1,s,
                               x,loc,si,ci)

        if s != -1.0:
            print 'collision', f,ci, s, x, loc
            x=numpy.array(x)

            c=gf.GetOutput().GetCell(ci)

            p0=numpy.array(c.GetPoints().GetPoint(0))
            p1=numpy.array(c.GetPoints().GetPoint(1))

            n=numpy.zeros(3)

            n[0]=(p1-p0)[1]
            n[1]=(p0-p1)[0]

            n=n/numpy.sqrt(sum(n**2))

            n=n*numpy.sign(numpy.dot(n,(p-pa)))

            print 'Before',  p0,p

            p=x+(1.0-s)*((p-pa)-(1.0+e)*n*(numpy.dot(n,(p-pa))))

            print 'After', pa,p,n

            theta=numpy.arccos(numpy.dot(n,(p-pa))/numpy.sqrt(numpy.dot(p-pa,p-p0)))

            coldat=collisionInfo(x,ci,theta,self.t+s*dt)

            if type(v) != type(None):
                v-=(1.0+e)*n*(numpy.dot(n,v))
                px,cr,vs=self.collide(p,v,x-1.0e-3*n,dt)
                p=px+x-1.0e-3*n
                return p-pa, coldat, vs-self.v
            else:
                p=x-1.0e-3*n+self.collide(p,None,x-1.0e-3*n,dt)[0]
                return p-pa, coldat
               
        if type(v) != type(None):
            return p-pa, None, v-self.v
        else:
            return p-pa, None
                               
        




reader.SetFileName('simple_0.vtu')
file1=reader.GetOutput()
file1.Update()

l1.SetDataSet(file1)
l2.SetDataSet(file2)
l1.BuildLocator()
l2.BuildLocator()


class particle_bucket(object):

    def __init__(self,X,V,t,dt,filename='data.dat'):

        self.particles=[]

        for x,v in zip(X,V):
            print x,v
            self.particles.append(particle(x,v,t,dt))
        self.t=t
        self.p=X
        self.v=V
        self.dt=dt
        self.file=open(filename,'w')

    def update(self):
        for p in self.particles:
            p.update()

        self.t+=self.dt

    def collisions(self):
        return itertools.chain(*[p.collisions for p in self.particles])


    def write(self):
        
        self.file.write('%f'%self.t)

        for p in self.p.ravel():
            self.file.write(' %f'%p)

        for v in self.v.ravel():
            self.file.write(' %f'%v)

        self.file.write('\n')

        

    

        



#X=numpy.array([[0.0,0.0,0.0],
#               [0.5,0.0,0.0],
#               [0.3,0.0,0.0]])
#V=numpy.array([[0.0,1.0,0.0],
#               [0.0,0.0,0.0],
#               [0.0,0.0,0.0]],)

N=1000

X=2.0*(numpy.random.random((N,3))-0.5)
X[:,2]=0
V=numpy.zeros((N,3))

pb=particle_bucket(X,V,0.0,0.01)
pb.write()


for k in range(300):
    print pb.t
    pb.update()
    pb.write()
pb.file.flush()




def plot(filename='data.dat'):
    f=open(filename,'r')

    t=[]
    x=[]
    y=[]
    u=[]
    v=[]

    for m in f.readlines():
        r=[float(M) for M in m.split()]


        print r
        t.append(r[0])
        n=len(r[1:])/6
        X=r[1::3]
        Y=r[2::3]
        print 
        x.append(X[:n])
        u.append(X[n:])
        y.append(Y[:n])
        v.append(Y[n:])

    
    print x
    print y

    x=numpy.array(x)
    y=numpy.array(y)

    p.plot(x,y)

        
