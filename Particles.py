import numpy
import TemporalCache
import IO
import DragModels
import Collision 

import vtk
import scipy.linalg as la

class particle(object):
    """Class representing a single Lagrangian particle with mass"""

    def __init__(self,p,v,t=0.0,dt=1.0,tc=None,u=numpy.zeros(3),
                 gp=numpy.zeros(3),rho=2.5e3,g=numpy.zeros(3),
                 omega=numpy.zeros(3),d=40e-6,boundary=None,
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
        self.boundary=boundary
        self.e=e
        self.drag=drag

    def update(self):
        """Update the state of the particle to the next time level

        The method uses relatively simple RK4 time integration."""

        collision=False

        kp1=self.v
        kv1=self.force(self.p,
                       self.v,
                       self.t)

        s,c1,v1=self.collide(kp1,0.5*self.dt,self.v,f=kv1)
        kp2=self.v+0.5*self.dt*kv1
        kv2=self.force(self.p+s,
                       self.v+v1,
                       self.t+0.5*self.dt)

        s,c2,v2=self.collide(kp2,0.5*self.dt,v=self.v,f=kv2)
        kp3=self.v+0.5*self.dt*kv2
        kv3=self.force(self.p+s,
                       self.v+v2,
                       self.t+0.5*self.dt)


        s,c3,v3=self.collide(kp3,self.dt,v=self.v,f=kv3)
        kp4=self.v+self.dt*kv3
        kv4=self.force(self.p+s,
                       self.v+v3,
                       self.t+self.dt)


        s,c4,v4=self.collide((kp1+2.0*(kp2+kp3)+kp4)/6.0,self.dt,
                             v=self.v,f=(kv1+2.0*(kv2+kv3)+kv4)/6.0)
        self.p+=s
        self.v+=v4
        if type(c4)!=type(None):
            self.collisions+=c4


        self.t+=self.dt
        self.u[:],self.gp[:]=self.picker(self.p,self.t)

    def force(self,p,v,t):
        
        """Calculate the sum of the forces on the particle.

        Args:
            p (float): Location at which forcing is evaluated.
            v (float): Particle velocity at which forcing is evaluated.
            t (float): Time at which particle is evaluated.
        """


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
            elif c.GetCellType()==vtk.VTK_TRIANGLE:
                ct=c
                ls=[0,1,2]
            elif c.GetCellType()==vtk.VTK_QUADRATIC_TETRA:
                ct=vtk.vtkTetra()
                ct.GetPoints().SetPoint(0,c.GetPoints().GetPoint(0))
                ct.GetPoints().SetPoint(1,c.GetPoints().GetPoint(1))
                ct.GetPoints().SetPoint(2,c.GetPoints().GetPoint(2))
                ct.GetPoints().SetPoint(3,c.GetPoints().GetPoint(3))
                ls=[0,1,2,3]
            elif c.GetCellType()==vtk.VTK_TETRA:
                ct=c
                ls=[0,1,2,3]

            N=ct.GetNumberOfPoints()-1
            x=numpy.zeros(ct.GetNumberOfPoints())
            args=[ct.GetPoints().GetPoint(i)[:N] for i in range(N+1)]
            args.append(x)
            ct.BarycentricCoords(p[:N],*args)

            collision=not Collision.testInCell(ct,p)          
        

            pids=c.GetPointIds()
            
            data_u=file.GetPointData().GetVectors('Velocity')
            data_p=file.GetPointData().GetScalars('Pressure')


            sf=numpy.zeros(c.GetNumberOfPoints())
            df=numpy.zeros(N*c.GetNumberOfPoints())
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


    def collide(self,k,dt,v=None,f=None,pa=None,level=0):

        """Collision detection routine.

        Args:
            k  (float): Displacement
            dt (float): Timestep
            v  (float, optional): velocity
            f  (float, optional): forcing
            pa (float, optional): starting position in subcycle
            level (int) count to control maximum depth
        """

        if type(pa)==type(None):
            pa=self.p

        if level==10 :
            if type(v) != type(None):
                return k*dt, None, v-self.v
            else:
                return k*dt, None
            


        p=pa+dt*k

        s=vtk.mutable(-1.0)
        x=[0.0,0.0,0.0]
        loc=[0.0,0.0,0.0]
        si=vtk.mutable(0)
        ci=vtk.mutable(0)

        intersect=self.boundary.bndl.IntersectWithLine(pa,p,
                               1.0e-8,s,
                               x,loc,si,ci)

        if s != -1.0:
            print 'collision', intersect,ci, s, x, loc
            x=numpy.array(x)

            c=self.boundary.bnd.GetCell(ci)

            p0=numpy.array(c.GetPoints().GetPoint(0))
            p1=numpy.array(c.GetPoints().GetPoint(1))

            n=numpy.zeros(3)

            n[0]=(p1-p0)[1]
            n[1]=(p0-p1)[0]

            n=n/numpy.sqrt(sum(n**2))

            n=n*numpy.sign(numpy.dot(n,(p-pa)))

            p=x+dt*(k-(1.0+self.e)*n*(numpy.dot(n,k)))


            theta=abs(numpy.arcsin(numpy.dot(n,(x-pa))/numpy.sqrt(numpy.dot(x-pa,x-pa))))

            coldat=[]

            


            if type(v) !=type(None):
                vs=v+s*dt*f
            else:
                vs=self.v+s*dt*f

            print 'Before',  p0,p,vs

            coldat.append(Collision.collisionInfo(x,vs,ci,theta,self.t+s*dt))
            vs+=-(1.0+self.e)*n*numpy.dot(n,vs)
            

            if type(v) != type(None):

                print 'After V1:', pa,p,n,vs,f


                px,cr,vo=self.collide(vs,(1-s)*dt,v=vs,f=f,pa=x+1.0e-10*vs,level=level+1)
                p=px+x+1.0e-10*vs

                print 'After V2:', pa,p,n,vs,f
                
                if cr:
                    coldat+=cr

                return p-pa, coldat, vo
            else:
                p=x+1.0e-8*vs+self.collide(vs,(1-s)*dt,f=f,pa=x+1.0e-8*vs,level=level+1)[0]


                print 'After', pa,p,n

                return p-pa, coldat
               
        if type(v) != type(None):
            
            return p-pa, None, v+dt*f-self.v
        else:
            return p-pa, None
                               

class particle_bucket(object):

    """Class for a container for multiple Lagrangian particles."""

    def __init__(self,X,V,t=0,dt=1.0e-3,filename=None,
                 base_name='',U=None,GP=None,rho=2.5e3,g=numpy.zeros(3),
                 omega=numpy.zeros(3),d=40.e-6,boundary=None,e=0.99,tc=None):
        """Initialize the bucket
        
        Args:
            X (float): Initial particle positions.
            V (float): Initial velocities
        """


        if tc:
            self.tc=tc
        else:
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
                                           d=d,boundary=boundary,e=e))
        self.t=t
        self.p=X
        self.v=V
        self.u=U
        self.gp=GP
        self.dt=dt
        self.boundary=boundary
        if filename: self.file=open(filename,'w')

    def update(self):
        """ Update all the particles in the bucket to the next time level."""
        self.tc.range(self.t,self.t+self.dt)
        for p in self.particles:
            p.update()

        self.t+=self.dt

    def collisions(self):
        """Collect all collisions felt by particles in the bucket"""
        return itertools.chain(*[p.collisions for p in self.particles])


    def write(self):
        """Write timelevel data to file."""
        
        self.file.write('%f'%self.t)

        for p in self.p.ravel():
            self.file.write(' %f'%p)

        for v in self.v.ravel():
            self.file.write(' %f'%v)

        self.file.write('\n')

    def run(self,t):
        """Drive particles forward until a given time."""
        while self.t<t:
            self.update()
            self.write()
        self.file.flush()
