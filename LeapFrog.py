"""
Module aims to allow LeapFrog integration on gravitational movement dependent primarily on distance

An acceleration function is necessary.
"""

def v0_half(v0, a, dt):
    # called to advance v0 half a step
    vh =  v0 + 0.5*a*dt
    return vh

def LeapFrog(r,v,a,dt, var):
    # to be called per timestep.
    # r,v,a different per call, 
    # var is the secondary parameter 
    # for a and can vary.
    R = r + v*dt
    V = v + a(R,var)*dt
    return R, V


