""" Module containing timestepping methods, all registered in the module dictionary variable 'methods' 

Methods have the generic signature 

    update_NAME(particle_model.Particle particle)
"""

from particle_model.Debug import profile, logger

def update_euler(self):
    """Update the state of the particle to the next time level
    
    The method uses the forward Euler method"""
    
    kap = (self.vel, self.force(self.pos,
                                self.vel,
                                self.time, drag=False), self.time)
    
    step, col, vel, col_vel = self.collide(self.vel,
                                           self.delta_t,
                                           self.vel,
                                           force=kap[1],
                                           drag=True)

    self.pos += step
    self.vel += vel

    if col:
        self.collisions += col

        kap = (self.vel, self.force(col[-1].pos+1.0e-8*col_vel,
                                    col_vel,
                                    col[-1].time+1.0e-8, drag=False), col[-1].time+1.0e-8)

    return kap, col

        
@profile
def update_ab2(self):
    """Update the state of the particle to the next time level
        
    The method uses the Adams Bashforth second order method"""

    if len(self._old)>=1:

        kap = (self.vel, self.force(self.pos,
                                    self.vel,
                                    self.time, drag=False), self.time)

        beta = 0.5*self.delta_t/(self.time-self.get_old(0,2))

        step, col, vel, col_vel = self.collide((1.0+beta)*self.vel-beta*self.get_old(0,0),
                                               self.delta_t,
                                               self.vel,
                                               force=(1.0+beta)*kap[1]-beta*self.get_old(0,1),
                                               drag=True)

        if col:
            self.collisions += col
            
            kap = (self.vel, self.force(col[-1].pos+1.0e-8*col_vel,
                                        col_vel,
                                        col[-1].time+1.0e-8, drag=False), col[-1].time+1.0e-8)

        self.pos += step
        self.vel += vel

        self.time += self.delta_t

    else:
        ## reduced to using the Euler method for the first timestep:

        kap, col = update_euler(self)

    self.set_old(kap,1)

    return kap, col

def update_ab3(self):
    """Update the state of the particle to the next time level

    The method uses the Adams Bashforth third order method"""

    if len(self._old)>=2:

        kap = (self.vel, self.force(self.pos,
                                    self.vel,
                                    self.time, drag=False), self.time)

        beta = -(1.0/6.0)*(self.delta_t*(5.0*self.delta_t+3.0*(self.time-self.get_old(0,2)))
                          /((self.time-self.get_old(0,2))*(self.get_old(0,2)-self.get_old(1,2))))
        gamma =  (1.0/6.0)*(self.delta_t*(2.0*self.delta_t+3.0*(self.time-self.get_old(0,2)))
                          /((self.time-self.get_old(1,2))*(self.get_old(0,2)-self.get_old(1,2))))

        step, col, vel, col_vel = self.collide((1.0-beta-gamma)*self.vel+beta*self.get_old(0,0)
                                               +gamma*self.get_old(1,0),
                                               self.delta_t,
                                               self.vel,
                                               force=(1.0-beta-gamma)*kap[1]+beta*self.get_old(0,1)
                                               +gamma*self.get_old(1,1),
                                               drag=True)


        if col:
            self.collisions += col

            kap = (self.vel, self.force(col[-1].pos+1.0e-8*col_vel,
                                        col_vel,
                                        col[-1].time+1.0e-8, drag=False), col[-1].time+1.0e-8)
            self.set_old(kap, 1)

        else:

            self.set_old(kap, 2)

        self.pos += step
        self.vel += vel

        self.time += self.delta_t

    else:
        ## reduced to using Adams Bashforth 2nd order method for the second timestep:

        tmp = self.get_old(0)
        kap, col = update_ab2(self)
        if tmp and not col:
            self._old = self._old + tmp

    return kap, col

def update_rk4(self):
    """Update the state of the particle to the next time level

    The method uses relatively simple RK4 time integration."""

    kap1 = (self.vel, self.force(self.pos,
                                 self.vel,
                                 self.time))

    step, col, vel, col_vel = self.collide(kap1[0], 0.5 * self.delta_t,
                                           self.vel,
                                           force=kap1[1])

    kap2 = (self.vel + 0.5*self.delta_t * kap1[1],
            self.force(self.pos + step,
                       self.vel + vel,
                       self.time + 0.5 * self.delta_t))

    step, col, vel, col_vel = self.collide(kap2[0], 0.5 * self.delta_t,
                                               vel=self.vel, force=kap2[1])
    kap3 = (self.vel+0.5 * self.delta_t * kap2[1],
            self.force(self.pos + step,
                       self.vel + vel,
                       self.time + 0.5*self.delta_t))

    step, col, vel, col_vel = self.collide(kap3[0], self.delta_t,
                                           vel=self.vel,
                                           force=kap3[1])
    kap4 = (self.vel + self.delta_t * kap3[1],
            self.force(self.pos + step,
                       self.vel + vel,
                       self.time + self.delta_t))

    step, col, vel, col_vel = self.collide((kap1[0]+2.0*(kap2[0]+kap3[0])+kap4[0])/6.0,
                                               self.delta_t, vel=self.vel,
                                               force=(kap1[1]+2.0*(kap2[1]+kap3[1])+kap4[1])/6.0)
    self.pos += step
    self.vel += vel
    if col:
        self.collisions += col


    self.time += self.delta_t

def update_rk2(self):
    """Update the state of the particle to the next time level

    The method uses second order Runge-Kutta time integration."""

    kap1 = (self.vel, self.force(self.pos,
                                 self.vel,
                                 self.time))

    step, col, vel, col_vel = self.collide(kap1[0], 0.5*self.delta_t,
                                           self.vel,
                                           force=kap1[1])

    kap2 = (self.vel + 0.5*self.delta_t * kap1[1],
            self.force(self.pos + 0.5*step,
                       self.vel + 0.5*vel,
                       self.time + 0.5 * self.delta_t))

    step, col, vel, col_vel = self.collide(kap2[0],
                                           self.delta_t, vel=self.vel,
                                           force=kap2[1])

    self.pos += step
    self.vel += vel
    if col:
        self.collisions += col


    self.time += self.delta_t

def update_rk3(self):
    """Update the state of the particle to the next time level

    The method uses third order Runge-Kutta time integration."""

    kap1 = (self.vel, self.force(self.pos,
                                 self.vel,
                                 self.time))

    step, col, vel, col_vel = self.collide(kap1[0], 0.5 * self.delta_t,
                                           self.vel,
                                           force=kap1[1])

    kap2 = (self.vel + 0.5*self.delta_t * kap1[1],
            self.force(self.pos + step,
                       self.vel + vel,
                       self.time + 0.5 * self.delta_t))

    step, col, vel, col_vel = self.collide(2.0*kap2[0]-kap1[0], self.delta_t,
                                               vel=self.vel, force=2.0*kap2[1]-kap1[1])
    kap3 = (self.vel+self.delta_t * (2.0*kap2[1]-kap1[1]),
            self.force(self.pos + step,
                       self.vel + vel,
                       self.time + self.delta_t))

    step, col, vel, col_vel = self.collide((kap1[0]+4.0*kap2[0]+kap3[0])/6.0,
                                           self.delta_t, vel=self.vel,
                                           force=(kap1[1]+4.0*kap2[1]+kap3[1])/6.0)

    self.pos += step
    self.vel += vel
    if col:
        self.collisions += col


    self.time += self.delta_t

methods = { "ForwardEuler":update_euler,
            "AdamsBashforth1":update_euler,
            "AdamsBashforth2":update_ab2,
            "AdamsBashforth3":update_ab3,
            "RungeKutta1":update_euler,
            "RungeKutta2":update_rk2,
            "RungeKutta3":update_rk3,
            "RungeKutta4":update_rk4}
