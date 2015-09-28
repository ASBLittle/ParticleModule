import numpy

class collisionException(Exception):
    pass

class collisionInfo(object):
    def __init__(self,x,v,cell,angle,t):
        self.x=x
        self.v=v
        self.cell=cell
        self.angle=angle
        self.t=t



def MclauryMassCoeff(ci,n=2):

    k=1.0
    H=1.0
    F_s=1.0
    F_B=1.0

    def f(t):
        if numpy.tan(t)>1.0/3.0:
            return numpy.cos(t)**2
        else:
            return numpy.sin(2.0*t)-3.0*numpy.sin(t)**2

    V=numpy.sqrt(numpy.sum(ci.v**2))

    return k*H*F_s*F_B*V**n*f(ci.angle)
    
    
