import numpy
import TemporalCache
import IO
import DragModels
import Collision 

import vtk
import scipy.linalg as la

class particle(object):

    def __init__(self,p,v,t=0.0,dt=1.0,tc=None,u=numpy.zeros(3),
                 gp=numpy.zeros(3),rho=2.5e3,g=numpy.zeros(3),
                 omega=numpy.zeros(3),d=40e-6,bndl=None,bnd=None,
                 e=0.99,drag=DragModels.turbulent_drag):
        
        self.p=p
        self.v=v
        self.t=t
        self.dt=dt
        self.collisions=[]
        self.tc=tc
        self.rho=2.0e3
        self.g=g
        self.omega=omega
        self.diameter=d
        self.u=u
        self.gp=gp
        self.bndl=bndl
        self.bnd=bnd
        self.e=e
        self.drag=drag

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

        drag=self.drag(u,v,self.diameter)

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

            if c.GetCellType()==vtk.VTK_QUADRATIC_TRIANGLE:
                ct=vtk.vtkTriangle()
                ct.GetPoints().SetPoint(0,c.GetPoints().GetPoint(0))
                ct.GetPoints().SetPoint(1,c.GetPoints().GetPoint(1))
                ct.GetPoints().SetPoint(2,c.GetPoints().GetPoint(2))
                ls=[0,1,2]
            else:
                ct=c
                ls=[0,1,2]

            collision=not ct.PointInTriangle(p,ct.GetPoints().GetPoint(0),
                              ct.GetPoints().GetPoint(1),
                              ct.GetPoints().GetPoint(2),1.0e-6)

            v1=[0.0,0.0]
            v2=[0.0,0.0]
            v3=[0.0,0.0]

            ct.ProjectTo2D(ct.GetPoints().GetPoint(0),
                          ct.GetPoints().GetPoint(1),
                          ct.GetPoints().GetPoint(2),
                          v1,v2,v3)

            x=[0.0,0.0,0.0]

            ct.BarycentricCoords(p[:2],
                                ct.GetPoints().GetPoint(0)[:2],
                                ct.GetPoints().GetPoint(1)[:2],
                                ct.GetPoints().GetPoint(2)[:2],
                                x)
        

            pids=c.GetPointIds()
            
            data_u=file.GetPointData().GetVectors('Velocity')
            data_p=file.GetPointData().GetScalars('Pressure')


            sf=numpy.zeros(c.GetNumberOfPoints())
            df=numpy.zeros(2.0*c.GetNumberOfPoints())
            c.InterpolateFunctions(x,sf)
            c.InterpolateDerivs(x,df)

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

        f=self.bndl.IntersectWithLine(pa,p,
                               1.0e-1,s,
                               x,loc,si,ci)

        if s != -1.0:
            print 'collision', f,ci, s, x, loc
            x=numpy.array(x)

            c=self.bnd.GetCell(ci)

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

            coldat=Collision.collisionInfo(x,v,ci,theta,self.t+s*dt)

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

    def __init__(self,X,V,t=0,dt=1.0e-3,filename=None,
                 base_name='',U=None,GP=None,rho=2.5e3,g=numpy.zeros(3),
                 omega=numpy.zeros(3),d=40.e-6,bndl=None,bnd=None,e=0.99):

        self.tc=TemporalCache.TemporalCache(base_name)
        self.particles=[]

        if U==None:
            U=[None for i in range(X.shape[0])]
        if GP==None:
            GP=[None for i in range(X.shape[0])]

        for x,v,u,gp in zip(X,V,U,GP):
            print x,v
            self.particles.append(particle(x,v,t,dt,tc=self.tc,u=u,gp=gp,
                                           rho=rho,g=g,omega=omega,
                                           d=d,bndl=bndl,bnd=bnd,e=e))
        self.t=t
        self.p=X
        self.v=V
        self.u=U
        self.gp=GP
        self.dt=dt
        if filename: self.file=open(filename,'w')

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

    def run(self,t):
        while self.t<t:
            self.update()
            self.write()
        self.file.flush()
