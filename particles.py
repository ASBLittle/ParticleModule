import vtk
import numpy
import scipy.linalg as la
import pylab as p
import itertools
import glob

class TemporalCache(object):
    
    def __init__(self,base_name):
        files=glob.glob(base_name+'*.vtu')

        self.data=[]

        self.lower=0
        self.upper=0

        

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
        
        self.open(0)
            
    def reset(self):
        self.lower=0
        self.upper=0

    def range(self,a,b):
        if self.data[self.lower][0]>a:
            self.reset()
        while self.lower<len(self.data)-2 and self.data[self.lower+1][0]<=a:
            self.close(self.lower)
            self.lower+=1
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
        if t1==t0: t1=1e300

        return self.data[lower:lower+2], (t-t0)/(t1-t0)

bnd_reader=vtk.vtkXMLUnstructuredGridReader()
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
    def __init__(self,x,v,cell,angle,t):
        self.x=x
        self.v=v
        self.cell=cell
        self.angle=angle
        self.t=t

class particle(object):

    def __init__(self,p,v,t=0.0,dt=1.0,tc=None,u=numpy.zeros(3),gp=numpy.zeros(3)):
        
        self.p=p
        self.v=v
        self.t=t
        self.dt=dt
        self.collisions=[]
        self.tc=tc
        self.rho=2.0e3
        self.g=numpy.array((0.,0.,0.))
        self.omega=numpy.array((0.,0.,100.0*numpy.pi))
        self.diameter=40.0e-6
        self.u=u
        self.gp=gp

    def update(self):
        # Simple RK4 integration

        collision=False

        kp1=self.v
        kv1=self.force(self.p,
                       self.v,
                       self.t)

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
        self.u[:],self.gp[:]=self.picker(self.p,self.t)

    def force(self,p,v,t):
        
        ## Drag term:


        u,grad_p=self.picker(p,t)

#        if collision:
#            raise collisionException


        beta=0.44*3.0/32.0/self.diameter*numpy.sqrt(sum((u-v)**2))

        drag=beta*(u-v)

        return grad_p/self.rho+drag+(-2.0*numpy.cross(self.omega,v)+self.g-numpy.cross(self.omega,numpy.cross(self.omega,p)))


    def picker(self,p,t):

        
        def fpick(file,locator):

            locator.BuildLocatorIfNeeded()


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
            
            data_u=file.GetPointData().GetVectors('Velocity')
            data_p=file.GetPointData().GetScalars('Pressure')


            sf=numpy.zeros(c.GetNumberOfPoints())
            c.InterpolateFunctions(p,sf)

            rhs=numpy.zeros(2)
            rhs[0]=data_p.GetValue(c.GetPointId(1))-data_p.GetValue(c.GetPointId(0))
            rhs[1]=data_p.GetValue(c.GetPointId(2))-data_p.GetValue(c.GetPointId(0))

            A=numpy.zeros((2,2))


            p0=numpy.array(c.GetPoints().GetPoint(0))
            p1=numpy.array(c.GetPoints().GetPoint(1))
            p2=numpy.array(c.GetPoints().GetPoint(2))
            A[0,:]=(p1-p0)[:2]
            A[1,:]=(p2-p0)[:2]

            A=la.inv(A)
            

            out=numpy.zeros(3)
            for k in range(c.GetNumberOfPoints()):
                out+=sf[k]*numpy.array(data_u.GetTuple(pids.GetId(k)))
                

            gp=numpy.zeros(3)
            gp[:2]=numpy.dot(A,rhs)


            return out,gp

        data,alpha=self.tc(t)

        u0,gp0=fpick(data[0][2],data[0][3])
        u1,gp1=fpick(data[1][2],data[1][3])
            
        return (1.0-alpha)*u0+alpha*u1,(1.0-alpha)*gp0+alpha*gp1


    def collide(self,p,v=None,pa=None,dt=None,level=0):
        if level==10 :
            if type(v) != type(None):
                return p-pa, None, v-self.v
            else:
                return p-pa, None
            
        if type(pa)==type(None):
            pa=self.p
            dt=self.dt
        e=0.95
        

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

            coldat=collisionInfo(x,v,ci,theta,self.t+s*dt)

            if type(v) != type(None):
                v-=(1.0+e)*n*(numpy.dot(n,v))
                px,cr,vs=self.collide(p,v,x-1.0e-3*n,dt,level=level+1)
                p=px+x-1.0e-3*n
                return p-pa, coldat, vs
            else:
                p=x-1.0e-3*n+self.collide(p,None,x-1.0e-3*n,dt,level=level+1)[0]
                return p-pa, coldat
               
        if type(v) != type(None):
            return p-pa, None, v-self.v
        else:
            return p-pa, None
                               

class particle_bucket(object):

    def __init__(self,X,V,t,dt,filename='data.dat',base_name='',U=None,GP=None):

        self.tc=TemporalCache(base_name)
        self.particles=[]


        for x,v,u,gp in zip(X,V,U,GP):
            print x,v
            self.particles.append(particle(x,v,t,dt,tc=self.tc,u=u,gp=gp))
        self.t=t
        self.p=X
        self.v=V
        self.u=U
        self.gp=GP
        self.dt=dt
        self.file=open(filename,'w')

    def update(self):
        self.tc.range(self.t,self.t+self.dt)
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

N=300

X=2.0*(numpy.random.random((N,3))-0.5)
X[:,2]=0
V=numpy.zeros((N,3))

pb=particle_bucket(X,V,0.0,0.01,base_name='simple')
pb.write()


for k in range(500):
    print pb.t
    pb.update()
    pb.write()
pb.file.flush()




def get_data(filename='data.dat'):
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

    f.close()

    return t,x,y

def plot(filename='data.dat'):

    t,x,y = get_data(filename)

    p.plot(x,y)

def movie(parbuck,filename='data.dat'):

    col=[c for c in pb.collisions()]

    col.sort(key=lambda x:x.t)

    t,x,y = get_data(filename)

    for k in range(len(t)):
        p.clf()
        p.plot(x[:k+1],y[:k+1],zorder=1)
        p.scatter(x[k],y[k],s=10,zorder=5)

        for c in col:
            if c.t<t[k]-0.01:
                p.scatter(c.x[0],c.x[1],c='r',s=30)
            elif c.t<=t[k]:
                p.scatter(c.x[0],c.x[1],c='k',s=50)
            else:
                break

        p.axis([0.05,0.15,-0.06,0.04])

        p.title('t=%4.3f'%t[k])

        p.savefig('movie/testF%05d.png'%k)
        
p.figure()

movie(pb)


