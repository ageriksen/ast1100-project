"""
Rewrite of first code to clarify
"""

##############################################
# imports
import numpy as np
import matplotlib.pyplot as plt
import MySolarSystem as M
import LeapFrog as lf
import scipy as scp


##############################################
# functions

def e_theta(theta, r):
    """
    transforms desired angle off of radial 
    vector to angle off of 1. axis.
    returns radians. arcatn2 input is (y, x)
    """
    theta = theta + np.arctan2(r[1], r[0])

def launchPosition(r, R, e_theta, theta):
    """
    position along surface of planet 0 to
    launch from.
    r is the position vector of the planet, 
    R is the radius of the planet, 
    e_theta is function above
    theta is angle from r in radians. 
    """
    return r + R*e_theta(theta, r)

def accelerate(r, m): 
    """
    calculates acceleration on body
    from gravitational pull from 
    body of mass m at a distance r
    """
    r_ = np.linalg.norm(r) # lenght of distance
    return -( G*m*r )/( r_**3 )

def planetvelocity(nr, time, epsilon): 
    """
    Approximating momentary velocity from 2 close
    positions around a time "time".
    nr is which planet,
    epsilon is the small variation across which we
    want to approximate velocity from. 
    """
    Dr = (
        pos_func(time+epsilon)[:,nr] 
        - pos_func(time-epsilon)[:,nr]
         ) # delta r
    Dt = 2.*epsilon # delta t
    return Dr / Dt

def journey(rs, vs, t, dt):
    """
    simulations of journey between planets. 
    rs is the satellite position vector. 
    vs is the satellite velocity vector. 
    t is the initial time of launch.
    dt is the distance in time per step.
    """
    # half-step of velocity:
    for j in range(N+1):
        # set v0 to v_{1/2}
        if j < N:
            vs[0] = v0_half( vs[0], accelerate( rs[0,:] - pos_func(t)[:,j] ), dt )
        else: 
            vs[0] = v0_half( vs[0], accelerate( rs[0,:], m[j] ), dt )
    
    # initiating variables for storing least distance and timestamp.
    rs1_m = 1 # AU
    ts1_m = 0 # yrs.
    # fulfilling the integration proper:
    i = 0
    maxDistance = np.linalg.norm(pos_func(t)[:,1]+r_stable+0.001)
    while np.linalg.norm(rs) < maxDistance:
    # iterate over planets for acceleration
        for j in range(N+1):
            #Leapfrog takes r, v, a, dt, var
            if j < N:
                # iterating over all the planets
                sat = (
                    LeapFrog( (rs[i] - pos_func(t)[:,j]), vs[i], accelerate,
                        dt, m[j]) )
                rs[i+1] += rs[i]+sat[0]
                vs[i+1] += vs[i]+sat[1]

            else:
                # pull of star.
                sat = LeapFrog( rs[i,:], vs[i], acelerrate, dt, m[j] )
                rs[i+1] += rs[i]+sat[0]
                vs[i+1] += vs[i]+sat[1]
        t += dt
        R = pod_func(t)[:,1] - rs[i+1]
        if linalg.norm(R) < linalg.norm(r_s1m):
            rs1_m = R
            ts1_m = t

###############################################################
# stating constants:
G = 4*np.pi**2 #AU**3yr**-2M_s**-1  gravitational constant
T = 10000 #1e4 K temperature in Kelvin
AU = 149597870700 #m in 1 AU

####################
# collectim simulation data from planet_orbits
T_max = 25 #yrs, timespan of simulation
n = 20000 # nr. of timesteps per year.
t_steps = n*T_max #timesteps for entire simulation

###################
#loading arrays
infile = open('positionsHomePlanet.npy', 'rb')
pos_p = load(infile)
time = np.linspace(0, T_max, t_steps) 

pos_func = scp.interpolate.interp1d(time, pos_p)
pos_p0 = pos_func(time[0])

#####################
# data from MySolarSystem
print '----------------------------------------------------------'
syst = M.Myseed() # instancing solar system with seed
N = M.Nplanets(syst) # nr. of planets in system
phi0 = M.init_angle(syst) # initial angles in shared orbital plane
m_s = M.starmass(syst) # mass of star [solar masses]
m_p = M.p_mass(syst) # mass of planet [solar masses]
m = np.append(m_p, m_s) # array of mass of sun and planets
R_p = M.p_radius(syst) #radius of planets[km]
R_p = R_p/AU #planet radius in AU

#########
#escape velocity
v_esc = sqrt( (2*G*m_p[0]) / R_p[0] ) 
