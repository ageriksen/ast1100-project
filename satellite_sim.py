"""
Rewrite of first code to clarify
table of contents:
header                approx line

imports                  20
functions                30
constants               140
planet_orbits           140
loading arrays          150
MySolarSystem           160
escape velocity         170
shortest distance       180
satellite launch        210
launch criteria         210
initial condition       210
launch                  220
plot                    240
"""

##############################################
# imports
import numpy as np
import numpy.linalg as nplin
import matplotlib.pyplot as plt
import MySolarSystem as M
import LeapFrog as lf
import scipy.interpolate as scp
import sys


##############################################
# functions

def e_theta(theta, r, i, j, c1=1, c2=1):
    """
    transforms desired angle off of radial 
    vector to angle off of 1. axis.
    returns radians. arcatn2 input is (y, x)
    """
    theta = theta + np.arctan2(r[1], r[0])
    e_t = np.zeros(2)
    e_t[0] = c1*i(theta)
    e_t[1] = c2*j(theta)
    return e_t

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

def journey(ms, rs, vs, t, dt, land='nan'):
    """
    simulations of journey between planets. 
    rs is the satellite position vector. 
    vs is the satellite velocity vector. 
    t is the initial time of launch.
    dt is the distance in time per step.
    """ 
    # half-step of velocity:
    a = np.zeros(2)
    r_p = pos_func(t)
    for j in range(N+1):
        if j < N:
            r = rs[0,:] - r_p[:,j]
            r_ = nplin.linalg.norm(r)
            a += - (G*m[j])/(r_**3) * r
        else:
            r_ = nplin.linalg.norm(rs[0])
            a += - (G*m[j])/(r_**3) * r
    vs[0] = vs[0] + 0.5*a*dt

    # initiating variables for storing least distance and timestamp.
    rs_min = np.zeros(2)
    t_min = 0 # yrs.
    i_min = 0
    first = 0
    # fulfilling the integration proper:
    for i in range(np.shape(rs)[0]-1):
    # iterate over planets for acceleration
        r_p = pos_func(t)
        a = np.zeros(2)
        rs[i+1] = rs[i] + vs[i]*dt
        for j in range(N+1):
            if j < N:
                r = rs[i] - r_p[:,j]
                r_ = nplin.linalg.norm(r)
                a += - (G*m[j])/(r_**3) * r
            else:
                r_ = nplin.linalg.norm(rs[i])
                a += - (G*m[j])/(r_**3) * r

        vs[i+1] = vs[i] + a*dt
        R = r_p[:,1] - rs[i+1]
        R_ = nplin.linalg.norm(R)
#        if abs(R_) < abs(nplin.linalg.norm(r_min)):
#            rs_min = R
#            t_min = t
#            i_min = i
        if land != 'nan':
            r_stable = nplin.linalg.norm(rs[i+1]) * np.sqrt( (m_s / m[1]) * (1./k) )
            if abs(R_) <= r_stable:
                if first < 1: 
                    print 'we can get a stable orbit'
                    i_min = i+1
                    rs_min = rs[i+1]
                    t_min = t
                    RS = nplin.linalg.norm(rs[i])
                    v_so = np.sqrt( (m[1]*G)/RS )
                    v_rel = v_so*e_theta(0, -R, np.sin, np.cos, -1)
                    v_p1 = planetvelocity(1, t, eps)
                    v_stable = v_p1 + v_rel
                    vs[i+1] = v_stable
                    first = 1
                else:
                    pass
                
        t += dt
    if land != 'nan':
        return rs, vs, t, rs_min, t_min, i_min
    else:
        return rs, vs, t

###############################################################
# stating constants:
G = 4*np.pi**2 #AU**3yr**-2M_s**-1  gravitational constant
T = 10000 #1e4 K temperature in Kelvin
AU = 149597870.7 #km in 1 AU
M_sol = 1.9891e30 #kg in 1 solar mass
m_sat = 1100 #kg placeholder satelite mass


#########
# collectim simulation data from planet_orbits
T_max = 25 #yrs, timespan of simulation
n = 20000 # nr. of timesteps per year.
t_steps = n*T_max #timesteps for entire simulation

#########
#loading arrays
path = '/home/anderger/Documents/Project_largeFiles/'
infile = open(path+'planetPositions.npy', 'rb')
pos_p, time = np.load(infile)
pos_func = scp.interpolate.interp1d(time, pos_p)

pos_p0 = pos_func(time[0])


#########
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
v_esc = np.sqrt( (2*G*m_p[0]) / nplin.linalg.norm(R_p[0]) ) 
print 'escaoe velocity: ', v_esc
print '----------------------------------------------------------'


################################################################
# estimating shortest distance between planets 0 and 1.
r_min = nplin.linalg.norm(abs(pos_p0[:,1] - pos_p0[:,0]))
r_ = nplin.linalg.norm(r_min)
t_min = 2 # initial value for least time, yrs
t = 3.7
dt = 1e-6
count = 0
#having narrowed search iteratively
while t < 3.71:
    r = pos_func(t)[:,1] - pos_func(t)[:,0]
    r_ = nplin.linalg.norm(r)
    if r_ < r_min:
        t_min = t
        r_min = r_
    t += dt
    count += 1

print '----------------------------------------------------------'
print 'least distance [AU]: ', r_min
print 'at time [yrs]: ', t_min

#approx stable orbit
k = 10
r_stable = pos_func(t_min)[:,1] * np.sqrt( (m_s / m[1]) * (1./k) )
print '----------------------------------------------------------'

#################################################################
#################################################################
# Satellite launch at planet 1

#########
# Launch criteria
theta = -np.pi*0.06925 # 45 degrees in radians
rp0 = pos_func(t_min)[:,0]
eps = 1e-6 # the small breadth epsilon for velocity
day = 1./365 # 1 day's portion of year 

# Initial conditions
#def e_theta(theta, r, i, j):
r0 = rp0 + R_p[0]*e_theta(np.pi/3, r, np.cos, np.sin)

v0 = planetvelocity( 0, t_min, eps ) + v_esc*e_theta( theta, rp0, np.cos, np.sin)*2


# Launch
#def journey(ms, rs, vs, t, dt):
#return rs, vs, t, rs1_m, ts1_m, i_min
length1 = int(np.shape(pos_p)[2]*0.5)

rs_launch = np.zeros( (length1, 2) )
vs_launch = np.zeros( (length1, 2) )
vs_launch[0] = v0
rs_launch[0] = r0

dt = 1e-9
t = t_min

print '----------------------------------------------------------'
print 'began launch at time t: ', t
#print 'initial acceleration'

rs_launch, vs_launch, t = journey( 
    m_sat, rs_launch, vs_launch, t, dt)
print 'finished launch at time t: ', t
print '----------------------------------------------------------'
#print '==============================='
#print 'vsl_launch', vs_launch
#print '==============================='

# mid flight
length2 = int(np.shape(pos_p)[2]*0.7765)
rs_mid = np.zeros( (length2, 2) )
vs_mid = np.zeros( (length2, 2) )
rs_mid[0] = rs_launch[-1]
vs_mid[0] = vs_launch[-1]

dt = 1e-7

print '----------------------------------------------------------'
print 'mid flight from time t: ', t
rs_mid, vs_mid, t = journey(
    m_sat, rs_mid, vs_mid, t, dt)
print 'to time t: ', t
print '----------------------------------------------------------'

#touchdown 
length3 = int(np.shape(pos_p)[2]*0.7)

rs_touch = np.zeros( (length3, 2) )
vs_touch = np.zeros( (length3, 2) )
vs_touch[0] = vs_mid[-1]
rs_touch[0] = rs_mid[-1]

dt = 1e-8

print '----------------------------------------------------------'
print 'began touch down at time t: ', t
#print 'initial acceleration'

rs_touch, vs_touch, t, rst_m, tst_m, ist_m = journey( 
    m_sat, rs_touch, vs_touch, t, dt, 1)
print 'finished touch down at time t: ', t
print '----------------------------------------------------------'

p0_start = pos_func(t_min)[:,0]
p1_last = pos_func(t)[:,1]

#####################################################################
# plotting journey
circ0 = plt.Circle( (p0_start[0],p0_start[1]), R_p[0], color='r')
circ1 = plt.Circle( (p1_last[0],p1_last[1]), R_p[1], color='r')
fig, ax = plt.subplots()
print '----------------------------------------------------------'
print 'visualizing'
#satelite
plt.plot( rs_launch[0,0], rs_launch[0,1], 'ro', label=('sat launch start') )
plt.plot( rs_launch[:,0], rs_launch[:,1], 'r' )
plt.plot( rs_mid[0,0], rs_mid[0,1], 'ro', label=('sat mid start') )
plt.plot( rs_mid[:,0], rs_mid[:,1], 'r' )
plt.plot( rs_touch[0,0], rs_touch[0,1], 'ro', label=('sat touchdown start') )
plt.plot( rs_touch[:,0], rs_touch[:,1], 'r' )
#least distance satellite
#plt.plot( rsl_m[0], rsl_m[1], color=('black'), marker=('v'), label=('least dist') )
plt.plot( pos_func(tst_m)[0,1], pos_func(tst_m)[1,1], 
    color=('black'), marker=('s'), label=('planet least') )
#radii of planets 0 and 1
ax.add_artist(circ0)
ax.add_artist(circ1)
#least distance planets
plt.plot( pos_func(t_min)[0,0], pos_func(t_min)[1,0], 'bo', 
    label=('planet 0 at time 3.71 yrs') )
plt.plot( pos_func(t_min)[0,1], pos_func(t_min)[1,1], 'bo', 
    label=('planet 1 at time 3.71 yrs') )
plt.plot( pos_func(t)[0,0], pos_func(t)[1,0], 'ro', 
    label=('planet 0 at time '+str(t)) )
plt.plot( pos_func(t)[0,1], pos_func(t)[1,1], 'ro', 
    label=('planet 1 at time '+str(t)) )
# planet orbits
plt.plot( pos_p[0, 0], pos_p[1, 0], 'g', label=('planet 0') )
plt.plot( pos_p[0, 1], pos_p[1, 1], 'g', label=('planet 1') )
for nr in range(2, N):
    plt.plot( pos_p[0, nr], pos_p[1, nr], 'y', label=('planet '+str(nr)) )
plt.legend()
plt.axis('equal')
plt.show()
print '----------------------------------------------------------'
