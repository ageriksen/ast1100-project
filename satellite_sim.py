"""
This program aims to plot the travel of a satelite from home planet
to target planet. 
"""

###############################################################################
# imports
from numpy import (
    array, pi, load, linspace, linalg,
    shape)
from matplotlib.pyplot import (
    plot, hold, show, legend)
import MySolarSystem as M
from LeapFrog import *
from scipy.interpolate import interp1d
import sys


###############################################################################
# functions

def accelerate(r, m):
    r_ = linalg.norm(r)
    return -(G*m*r)/(r_**3)

def velocity(dp, dti, scale='m/s'):
    'returns velocity in AU/yr or m/s'
    if scale == 'AU/yr':
        return dp/dt
    elif scale == 'm/s':
        dp = dp*AU
        dt = dp*yrs
        return dp/dt
    else:
        print 'wrong scale, buddy.'
        sys.exit

def e_theta(theta):
    return array((cos(theta), sin(theta)))

def launchPosition(r, R, e_theta, theta):
    """
    r[0] is the position vector of planet 0
    R[0] is the radius of planet 0.
    e_theta is a vector describing the 
    angle of the launch relative to r[0]
    """
    launch_pos = r + R*e_theta(theta)

###############################################################################
# set constants
G = 4*pi**2 # grav. constant
T = 10000# 1e4K temperature
#k # boltzmann constant

###############################################################################
# collecting from my system
syst = M.Myseed() # instancing SolarSystem with my seed.
N = M.Nplanets(syst) # nr. of planets
r0 = [[syst.x0[i],syst.y0[i]] for i in range(N) ]
v0 = [[syst.vx0[i],syst.vy0[i]] for i in range(N) ]
phi0 = M.init_angle(syst) # initial angle in shared orbital plane
m_s = M.starmass(syst) #starmass
m_p = M.p_mass(syst) # planetmass


###############################################################################
# collecting sim-data from planet_orbits
T_max = 25 #yrs. span of sim.
n = 20000 # nr. of timesteps in a year.
t_steps = n*T_max # timesteps for sim.
inFile = open('positionsHomePlanet.npy', 'rb')
pos_p = load(inFile)
time = linspace(0, T_max, t_steps) #bruk tall fra part2 t_steps, t_max

print shape(pos_p)

pos_func = interp1d(time, pos_p)
pos_p0 = pos_func(time[0])

###############################################################################
# figure time of shortest distance between the planets
r_min = linalg.norm(pos_p0[:,1])
for t in time:
    r = pos_func(t)[:,1] - pos_func(t)[:,0]
    r_ = linalg.norm(r)
    if r_ < r_min:
        t_min = t
        r_min = r_

#print 'least distance:'
#print r_min
#print 'at time: '
#print t_min




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

###############################################################################
# launch satellite slightly in front of planet 1, redo 

# half-step
#for i in planets:
#    for j in planets:
#        if i != j:
#            r[0] = r[0] + v0_half(
#                v[0,i,:], accelerate( (r[0,i,;]-r[0,j,:]), m), dt)
#            
## integrating
#for t in time:
#    for i in planets:
#       for j in planets:
#           if i != j:
#               r[t] = r[t] + ieapFrog(
#                    r[t,i,:]-r[t,j,:], v[t,i,:], accelerate, dt, m) 
#
###############################################################################
# untill satisfied/boost for proper velocity and direction. 


