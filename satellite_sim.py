"""
This program aims to plot the travel of a satelite from home planet
to target planet. 
"""

# imports
from numpy import (
    array, pi, load, linspace, linalg,
    shape)
from matplotlib.pyplot import (
    plot, hold, show)
import MySolarSystem as M
from LeapFrog import *
from scipy.interpolate import interp1d

# set constants
G = 4*pi**2 # grav. constant
T = 10000# 1e4K temperature
#k # boltzmann constant

# retrieving heavy planets
syst = M.Myseed() # instancing SolarSystem with my seed.
N = M.Nplanets(syst) # nr. of planets
r0 = [[syst.x0[i],syst.y0[i]] for i in range(N) ]
v0 = [[syst.vx0[i],syst.vy0[i]] for i in range(N) ]
phi0 = M.init_angle(syst) # initial angle in shared orbital plane
m_s = M.starmass(syst) #starmass
m_p = M.p_mass(syst) # planetmass


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



print 'visualising'
hold('on')
for nr  in range(N):
    plot(pos_p[0,nr], pos_p[1,nr], label=('planet ' +str(nr)))
#    plot(pos_func(t_min)[:,nr], label=('planet '+str(nr) + 'at time '+str(t_min)+'yrs'))
show()

# launch satellite slightly in front of planet 1, redo 
# untill satisfied/boost for proper velocity and direction. 


