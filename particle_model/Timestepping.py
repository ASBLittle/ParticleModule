""" Module containing timestepping methods, all registered in the module
dictionary variable 'methods'.

Methods have the generic signature

    update_NAME(particle_model.Particle particle)
"""

from particle_model.Debug import profile, logger
from particle_model import Collision

def update_euler(self, delta_t=None):
    """Update the state of the particle to the next time level

    The method uses the forward Euler method"""

    delta_t = delta_t or self.delta_t

    kap = (self.vel, self.force(self.pos,
                                self.vel,
                                self.time, drag=False), self.time)

    pos = self.pos+delta_t*self.vel
    vel = self.vel+delta_t*kap[1]

    try:
        self.pos, self.vel = self.check_collision_full(pos, self.pos,
                                                       vel, self.vel,
                                                       delta_t, drag=True)
    except Collision.CollisionException as col:
        vel = self.vel+col.delta_t*kap[1]
        C, fvel = self.drag_coefficient(col.pos, vel, self.time+col.delta_t, nearest = True)
        col.vel = (self.vel+col.delta_t*(kap[1]+C*fvel))/(1.0+col.delta_t*C)
        raise col
        
    self.time += delta_t

    return kap

     
@profile
def update_ab2(self, delta_t=None):
    """Update the state of the particle to the next time level

    The method uses the Adams Bashforth second order method"""

    delta_t = delta_t or self.delta_t

    if len(self._old) >= 1:

        kap = (self.vel, self.force(self.pos,
                                    self.vel,
                                    self.time, drag=False), self.time)

        beta = 0.5*self.delta_t/(self.time-self.get_old(0, 2))
        
        pos = self.pos+delta_t*((1+beta)*self.vel-beta*self.get_old(0, 0))
        vel = self.vel+delta_t*((1+beta)*kap[1]-beta*self.get_old(0, 1))

        try:
            self.pos, self.vel = self.check_collision_full(pos, self.pos,
                                                           vel, self.vel,
                                                           delta_t, drag=True)
        except Collision.CollisionException as col:
            beta = 0.5*col.delta_t/(self.time-self.get_old(0, 2))
            vel = self.vel+col.delta_t*(1+beta)*kap[1]-beta*self.get_old(0, 1)
            C, fvel = self.drag_coefficient(col.pos, vel, self.time+col.delta_t, nearest = True)
            col.vel = (self.vel+col.delta_t*(kap[1]+C*fvel))/(1.0+col.delta_t*C)
            raise col
        
        self.time += delta_t

    else:
        ## reduced to using the Euler method for the first timestep:

        kap = update_euler(self)

    self.set_old(kap, 1)

    return kap

def update_ab3(self, delta_t=None):
    """Update the state of the particle to the next time level

    The method uses the Adams Bashforth third order method"""

    delta_t = delta_t or self.delta_t

    if len(self._old) >= 2:

        kap = (self.vel, self.force(self.pos,
                                    self.vel,
                                    self.time, drag=False), self.time)

        beta = -(1.0/6.0)*(self.delta_t*(5.0*self.delta_t+3.0*(self.time-self.get_old(0, 2)))
                           /((self.time-self.get_old(0, 2))*(self.get_old(0, 2)-self.get_old(1, 2))))
        gamma =  (1.0/6.0)*(self.delta_t*(2.0*self.delta_t+3.0*(self.time-self.get_old(0, 2)))
                            /((self.time-self.get_old(1, 2))*(self.get_old(0, 2)-self.get_old(1, 2))))

        pos = self.pos+(1-beta-gamma)*self.vel+beta*self.get_old(0, 0)+gamma*self.get_old(1, 0)
        vel = self.vel+(1-beta-gamma)*kap[1]+beta*self.get_old(0, 1)+gamma*self.get_old(1, 0)

        try:
            self.pos, self.vel = self.check_collision_full(pos, self.pos,
                                                           vel, self.vel,
                                                           delta_t, drag=True)
        except Collision.CollisionException as col:
            beta = -(1.0/6.0)*(col.delta_t*(5.0*col.delta_t+3.0*(self.time-self.get_old(0, 2)))
                               /((self.time-self.get_old(0, 2))*(self.get_old(0, 2)-self.get_old(1, 2))))
            gamma =  (1.0/6.0)*(col.delta_t*(2.0*col.delta_t+3.0*(self.time-self.get_old(0, 2)))
                                /((self.time-self.get_old(1, 2))*(self.get_old(0, 2)-self.get_old(1, 2))))
            vel = self.vel+(1-beta-gamma)*kap[1]+beta*self.get_old(0, 1)+gamma*self.get_old(0, 1)
            C, fvel = self.drag_coefficient(col.pos, vel, self.time+col.delta_t, nearest = True)
            col.vel = (self.vel+col.delta_t*(kap[1]+C*fvel))/(1.0+col.delta_t*C)
            raise col

        self.set_old(kap, 2)

        self.time += self.delta_t

    else:
        ## reduced to using Adams Bashforth 2nd order method for the second timestep:

        try:
            tmp = [self.get_old(0)]
        except IndexError:
            tmp = []
        kap = update_ab2(self)
        if tmp:
            self._old = self._old + tmp

    return kap

def update_rk4(self, delta_t=None):
    """Update the state of the particle to the next time level

    The method uses relatively simple RK4 time integration."""

    delta_t = delta_t or self.delta_t

    try:
    
        kap1 = (self.vel, self.force(self.pos,
                                     self.vel,
                                     self.time))

        pos = self.pos+0.5*delta_t*kap1[0]
        vel = self.vel+0.5*delta_t*kap1[1]
        self.check_collision_full(pos, self.pos,
                                  vel, self.vel,
                                  0.5*delta_t, drag=False)
    
        kap2 = (self.vel + 0.5*delta_t*kap1[1],
                self.force(pos, vel, self.time + 0.5*delta_t))

        pos = self.pos+0.5*delta_t*kap2[0]
        vel = self.vel+0.5*delta_t*kap2[1]
        self.check_collision_full(pos, self.pos,
                                  vel, self.vel,
                                  0.5*delta_t, drag=False)


        kap3 = (self.vel+0.5*delta_t*kap2[1],
                self.force(pos, vel, self.time+0.5*delta_t))

        pos = self.pos+0.5*delta_t*kap3[0]
        vel = self.vel+0.5*delta_t*kap3[1]
        self.check_collision_full(pos, self.pos,
                                  vel, self.vel,
                                  0.5*delta_t, drag=False)


        kap4 = (self.vel + delta_t * kap3[1],
                self.force(pos, vel, self.time + delta_t))

        pos = self.pos+delta_t*(kap1[0]+2.0*kap2[0]+2.0*kap3[0]+kap4[0])/6.0
        vel = self.vel+delta_t*(kap1[1]+2.0*kap2[1]+2.0*kap3[1]+kap4[1])/6.0
        self.check_collision_full(pos, self.pos,
                                  vel, self.vel,
                                  delta_t, drag=False)

        self.pos = pos
        self.vel = vel

    except Collision.CollisionException as col:
        col.vel = self.vel+col.delta_t*kap1[0]
        raise col

    self.time += delta_t

def update_rk2(self, delta_t=None):
    """Update the state of the particle to the next time level

    The method uses second order Runge-Kutta time integration."""

    delta_t = delta_t or self.delta_t

    try:

        kap1 = (self.vel, self.force(self.pos,
                                     self.vel,
                                     self.time))

        pos = self.pos+0.5*delta_t*kap1[0]
        vel = self.vel+0.5*delta_t*kap1[1]
        self.check_collision_full(pos, self.pos,
                                  vel, self.vel,
                                  0.5*delta_t, drag=False)

        kap2 = (vel, self.force(pos, vel, self.time+0.5*delta_t))

        pos = self.pos+delta_t*kap2[0]
        vel = self.vel+delta_t*kap2[1]
        self.check_collision_full(pos, self.pos,
                                  vel, self.vel,
                                  delta_t, drag=False)

        self.pos = pos
        self.vel = vel

    except Collision.CollisionException as col:
        col.vel = self.vel+col.delta_t*kap1[0]
        raise col

    self.time += self.delta_t

def update_rk3(self, delta_t=None):
    """Update the state of the particle to the next time level

    The method uses third order Runge-Kutta time integration."""

    delta_t = delta_t or self.delta_t

    try:

        kap1 = (self.vel, self.force(self.pos,
                                     self.vel,
                                     self.time))

        pos = self.pos+0.5*delta_t*kap1[0]
        vel = self.vel+0.5*delta_t*kap1[1]

        self.check_collision_full(pos, self.pos,
                                  vel, self.vel,
                                  0.5*delta_t, drag=False)

        kap2 = (vel, self.force(pos, vel, self.time+0.5*delta_t))

        pos = self.pos+delta_t*(2.0*kap2[0]-kap1[0])
        vel = self.vel+delta_t*(2.0*kap2[1]-kap1[1])

        self.check_collision_full(pos, self.pos,
                                  vel, self.vel,
                                  delta_t, drag=False)

        kap3 = (vel, self.force(pos, vel, self.time+delta_t))

        pos = self.pos+delta_t*(kap1[0]+4.0*kap2[0]+kap3[0])/6.0
        vel = self.vel+delta_t*(kap1[1]+4.0*kap2[1]+kap3[1])/6.0

        self.pos, self.vel = self.check_collision_full(pos, self.pos,
                                                       vel, self.vel,
                                                       delta_t, drag=False)
    
    except Collision.CollisionException as col:
        col.vel = self.vel+col.delta_t*kap1[0]
        raise col

    self.time += self.delta_t

methods = {"ForwardEuler":update_euler,
           "AdamsBashforth1":update_euler,
           "AdamsBashforth2":update_ab2,
           "AdamsBashforth3":update_ab3,
           "RungeKutta1":update_euler,
           "RungeKutta2":update_rk2,
           "RungeKutta3":update_rk3,
           "RungeKutta4":update_rk4}
