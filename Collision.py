class collisionException(Exception):
    pass

class collisionInfo(object):
    def __init__(self,x,v,cell,angle,t):
        self.x=x
        self.v=v
        self.cell=cell
        self.angle=angle
        self.t=t
