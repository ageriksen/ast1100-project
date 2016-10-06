"""
This program aims to plot the travel of a satelite from home planet
to target planet. 
"""

###############################################################################
# imports
from numpy import (
    array, pi, cos, sin, sqrt,  load, linspace, zeros, linalg,
    shape, append)
from matplotlib.pyplot import (
    plot, hold, show, legend)
import MySolarSystem as M
from LeapFrog import *
from scipy.interpolate import interp1d
import sys


###############################################################################
# functions

def accelerate(r, m):
    """
    Function takes in the distance and mass of a body and
'   calculates the resulting acceleration another body would have towards it. 
    """
    r_ = linalg.norm(r)
    return -(G*m*r)/(r_**3)


def e_theta(theta):
    """
    function takes an angle theta and returns the cartesian
    unitvector from the polar angle. 
    """
    return array((cos(theta), sin(theta)))

def launchPosition(r, R, e_theta, theta):
    """
    r[0] is the position vector of planet 0
    R[0] is the radius of planet 0.
    e_theta is a vector describing the 
    angle of the launch relative to r[0]
    """
    return r + R*e_theta(theta)

def planetvelocity(nr, time, epsilon):
    """
    Aims to approximate momentary velocity from 2 close positions in time. 
    nr = which planet we are looking at
    time is the time we wish to find the velocity at
    epsilon is the small time variation we look at to estimate the velocity.
    """
    Dr = pos_func(time+epsilon)[:,nr] - pos_func(time-epsilon)[:,nr]
    Dt = 2*epsilon
    return Dr / Dt

def journey(rs, vs, t, dt):
    """
    Cronichles the journey of the sattelite through the solar system, 
    the starttime t and the length of each timestep.
    """
    for j in range(N+1):
        if j < N: 
            vs[0] = v0_half(vs[0], accelerate(rs[0,:]-pos_func(t)[:,j], m[j]), dt)
        else:
            vs[0] = v0_half(vs[0], accelerate(rs[0,:], m[j]), dt)
    
    for i in range(shape(rs)[1]):
        for j in range(N+1): #iterate over planets and star
            if j < N:
                #iterating over planets
                # LeapFrog(r, v, a, dt, var)
                rs[i+1], vs[i+1] = (
                    rs[i] + LeapFrog( (rs[i] - pos_func(t)[:,j]), vs[i], accelerate, 
                        dt, m[j])[0],
                    vs[i] + LeapFrog( (rs[i] - pos_func(t)[:,j]), vs[i], accelerate, 
                        dt, m[j])[1]
                                   )
            else:
                #pull from star
                rs[i+1], vs[i+1] = (
                    rs[i] + LeapFrog( rs[i,:], vs[i], accelerate, 
                        dt, m[j])[0],
                    vs[i] + LeapFrog( rs[i,:], vs[i], accelerate, 
                        dt, m[j])[1]
                                  )
        t += dt
        if linalg.norm(rs[i+1]) > 1:
            sys.exit 
    return rs, vs, t

###############################################################################
# set constants
G = 4*pi**2 # grav. constant
T = 10000# 1e4K temperature
AU = 149597870700 #m in 1 AU
#k # boltzmann constant

###############################################################################
# collecting from my system
print '---------------------------------------------------------------'
syst = M.Myseed() # instancing SolarSystem with my seed.
N = M.Nplanets(syst) # nr. of planets
phi0 = M.init_angle(syst) # initial angle in shared orbital plane
m_s = M.starmass(syst) #starmass
m_p = M.p_mass(syst) # planetmass
m = append(m_p, m_s) # new array of sun and planets
R_p = M.p_radius(syst) # planet radius
R_p = R_p/AU #planet radius in AU

v_esc = sqrt((2*G*m_p[0])/R_p[0]) # escape velocity og planet 0

###############################################################################
# collecting sim-data from planet_orbits
T_max = 25 #yrs. span of sim.
n = 20000 # nr. of timesteps in a year.
t_steps = n*T_max # timesteps for sim.
inFile = open('positionsHomePlanet.npy', 'rb')
pos_p = load(inFile)
time = linspace(0, T_max, t_steps) #bruk tall fra part2 t_steps, t_max

#print 'pos_p[2]'
#print shape(pos_p)[2]
#print type(shape(pos_p)[2])

pos_func = interp1d(time, pos_p)
pos_p0 = pos_func(time[0])

print '---------------------------------------------------------------'
###############################################################################
# figure time of shortest distance between the planets
r_min = linalg.norm(pos_p0[:,1])
r_ = linalg.norm(r_min)
for t in time:
    if t < 3:
        continue
    elif t > 4:
        break
    else:
        r = pos_func(t)[:,1] - pos_func(t)[:,0]
        r_ = linalg.norm(r)
        if r_ < r_min:
            t_min = t
            r_min = r_
            
"""consider saving this file as .npy and retrieveing it to save time, or
restrict the range of the loop to ~3yrs, or making it a while loop for 
easier manipulation"""
print '---------------------------------------------------------------'
print 'least distance [AU]:'
print r_min
print 'at time [yrs]: '
print t_min
print '---------------------------------------------------------------'







###############################################################################
# launch satellite slightly in front of planet 1, redo 

#launch criteria:
theta = 30 #degrees
rp0 = pos_func(t_min)[:,0]
epsilon = 1e-5

r0 = launchPosition(rp0, linalg.norm(pos_func(t_min)), e_theta, theta)
v0 = 2*v_esc*e_theta(theta) + planetvelocity(0, t_min, epsilon)
#v0 = planetvelocity(1, t_min, epsillon) - planetvelocity(0,t_min, epsilon)

# position for transit
length1 = int(shape(pos_p)[2]*0.5)

rs_launch = zeros((length1, 2))
vs_launch = zeros((length1, 2))
vs_launch[0] = v0
rs_launch[0] = r0

#dt = somenumber
#t = t_min
#rs, vs = journey(rs, vs, t, dt)
t = t_min
dt_m = 1e-6 #yrs timestep launch / landing
dt_l = 1e-3 #yrs timestep free flight


rs_launch, vs_launch, t = journey(rs_launch, vs_launch, t, dt_m)


####
#print 'visualising'
#plot(pos_func(t_min)[0,0], pos_func(t_min)[1,0],
#   'o', label=('planet 0 at time 3.71 yrs'))
#plot(pos_func(t_min)[0,1], pos_func(t_min)[1,1],
#   'o', label=('planet 1 at time 3.71 yrs'))
#hold('on')
#for nr  in range(N):
#    plot(pos_p[0,nr], pos_p[1,nr], label=('planet ' +str(nr)))
#legend()
#show()
####

###############################################################################
# untill satisfied/boost for proper velocity and direction. 
###############################################################################
"""
--------------------------------------------------------------
Run example:
>>> python satellite_sim.py 
---------------------------------------------------------------
myseed 80101 <type 'int'>
homeplanet is planet 0
number of planets 8
no index given
initial angle of planets' orbits, rad:
[ 0.          3.18581766  3.54960507  2.68861694  1.07636256  0.9163174
  4.38663014  2.8985559 ]
mass of star, [solar masses]
mass of planets, [solar masses]
no index given
radi of planets, [km]
---------------------------------------------------------------
---------------------------------------------------------------
least distance [AU]:
0.30641235698
at time [yrs]: 
3.70855741711
---------------------------------------------------------------
"""
