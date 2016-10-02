"""
exercise 1A.7:
builds on 1A.6. 
find temperature of gass and density to reach
escape vel. in 20 mins with fuel. 
can use to get full vid of launch in MCAst to check if successful. 

For now, copying major parts of 1A6
"""
#################################################################################
#imports
from AST1100SolarSystem import AST1100SolarSystem
from numpy import(
pi, sqrt, random, zeros, linspace, fft, histogram, mean, sum, where, nonzero, shape
)
from matplotlib.pyplot import plot, show, figure, bar, title, savefig
from MySolarSystem import * # module for all SolarSystem functions.
import time

#       #noting processor clock 
start_time = time.clock()

#       #retrieving solar system data
seed = Myseed() # fetching my seed from textfile myseed.txt through module. 
m_0 = p_mass(seed, 0) #Solar Masses, mass of home planet (home is 0)
R = p_radius(seed, 0) #km, radius of home planet (home is 0)
R = R*1000 # m, radius of planet.

        #setting important constants
AU = 149597870700 #m exactly defined [2]
G = 6.67408e-11 # N*m^2*kg^-2 Newton's gravitational constant [1]
mS = 1.988500e30 #kg solar mass[3]
m_p = m_0*mS # mass of planet 0 in kg
m_sat = 1000 #kg mass of sattelite, fuel mass ignored
#Boltzmann constant
k = 1.38064852e-23 # J*K^-1 Boltzamann constant [1]


        #motor:cubic box 
        #lenght:
L = 0.000001#1e-6 #m
#temperature inside
T = 10000 #1e4 #K
#nr. of particles
N = 100000 #1e5 #dimless
#mass per particle
m_h = 3.3476467e-27 # kg mass of Deuterium, aka. hydrogen gas or H_2


        #### 1: calculate escape velocity v_esc ####

# escape velocity is the minimal velocity to escape the gravitational pull
# of the homeplanet 0. this is when Ek > E_grav.
v_esc = sqrt((2*G*m_p)/R) # m/s
print '---------------------------------------------------------'
print 'escape velocity is:', v_esc, 'm/s'
print '---------------------------------------------------------'

#        #### 2: dedicate vectors for particles in box ####
#random.seed(100) # seeding the generator to retain numbers between runs
#l = 0
#h = L
#siz = N
#r = random.uniform(l,h,(siz,3))
#
#mu = 0
#o = sqrt((k*T)/m_h) # sigma, the standard deviation
#v = random.normal(mu,o,(siz,3))
#
#
#        #### test to check uniform and normal distribution ###
#        #hist, bins = histogram(r, 'auto')
#        #width = 0.7*(bins[1] - bins[0])
#        #center = (bins[:-1]+bins[1:])/2
#        #bar(center, hist, align='center', width=width)
#        #title('position')
#        #savefig('testfigures/histogram_test_position.pdf')
#        #show()
#        #hist, bins = histogram(v, 'auto')
#        #width = 0.7*(bins[1] - bins[0])
#        #center = (bins[:-1]+bins[1:])/2
#        #bar(center, hist, align='center', width=width)
#        #title('velocity')
#        #savefig('testfigures/histogram_test_velocity.pdf')
#        #show()
#        #### Test successfull. r is uniform and v is gaussian ###
#
#
#print 'testing generated values'
#print '---------------------------------------------------------'
#
#        #### 3a: is Ek = (3/2)*k*T? ####
#
#Ek_p = (3./2)*N*k*T
#Ek = (0.5*m_h)*sum( v**2 )
#rel_err = abs( (Ek - Ek_p)/Ek_p )
#print 'kinetic energy'
#print 'numerical', Ek_p
#print 'analytical', Ek
#print 'rel_err', rel_err
#print '---------------------------------------------------------'
#
#
#        #### 3b:test that the mean abs velocity follows it's relation ####
#
#absv_b = 2*N*sqrt( (2*k*T)/(pi*m_h) ) # <v> = integ (v*P(v)) | scaled up to attempt to reduce roundoff
#absv = sum(sqrt( sum( v**2, axis=1 ) )) #total absolute velocity
#relerr_absv = abs((absv- absv_b) / absv_b)
#print 'velocity'
#print 'analytical', absv_b
#print 'numerical', absv
#print 'rel_error', relerr_absv
#print '---------------------------------------------------------'
#
#
#        ### 4: simulate particle motion over ns with elastic collision of walls ###
#print '---------------------------------------------------------'
#print 'simulating particle motion in time'
#print '---------------------------------------------------------'
#Dt = 1e-9 #s, 1ns total time elapsed
#n = 1000 # number of timesteps within 1ns.
#def rocket(n, Dt,r,v, L):
#    dt = Dt/n # timestep
#    p_out = 0
#    n_out = 0
#    p_z = 0
#    for timestep in range(n):
#        for xyz in range(3):
#            roof = where(r[:,xyz] > L, 1,0)
#            floor = where(r[:,xyz] < 0, 1,0)
#            crash = roof+floor
#            if xyz == 2:
#                # placing hole in floor of z
#                n_x = where(
#    	    (r[nonzero(floor),0]>L/4.) & (r[nonzero(floor),0]<3.*L/4), 1,0)
#                n_y = where(
#    	    (r[nonzero(floor),1]>L/4.) & (r[nonzero(floor),1]<3.*L/4), 1,0)
#                n_pout = n_x*n_y
#                n_pin = where(n_x*n_y > 0, 0,1)
#                p_out += 2.*m_h*sum(abs(v[nonzero(n_pout),2]))
#                n_out += sum(n_pout)
#                r[nonzero(n_pout), xyz] += L
#                r[nonzero(n_pin),xyz] = - r[nonzero(n_pin),xyz]
#                r[nonzero(roof), xyz] = 2*L - r[nonzero(roof),xyz]
#                v[nonzero(n_pin),xyz] = -v[nonzero(n_pin),xyz]
#            else:
#                v[nonzero(crash),xyz] = -v[nonzero(crash),xyz]
#                r[nonzero(roof), xyz] = 2*L - r[nonzero(roof),xyz]
#                r[nonzero(floor), xyz] = -r[nonzero(floor), xyz]
#        p_z += 2.*m_h*sum(abs(v[nonzero(roof),2]))
#        r += v*dt
#    return p_z, r, v
#print 'motion simulation'
#print '---------------------------------------------------------'
#
#
#############################################################################################
####################                                         	############################
####################	THIS IS WHERE EXERCISE 1A.7 STARTS!!	############################
####################                                         	############################
#############################################################################################
#print '---------------------------------------------------------'
#print '---------------------------------------------------------'
#print 'exercise 1A.7'
#print '---------------------------------------------------------'
#print '---------------------------------------------------------'
#	### 1) calculate acceleration from eq. of state?
## eq. of stat: P = nkT | P=pressure; n=particle density; k=Boltzsmanns constant; T=temp
## numerical pressure: P = p/(A*Dt) | P=pressure; p=momentum; A=area; Dt=Delta t->timespan.
#
## could write: P = nkT = NkT/V = p/(A*Dt)
## | N=nr of particles; k=Boltzmann, T=temp; V=volume; p=momentum; A=area, Dt=delta t
## Assuming uniform pressure. Assuming constant pressure. with random particle position and
## direction of velocity, equally likely to hit a point in the hole as outside. 
## p_floor/(A_floor*Dt) = p_hole / (A_hole*Dt)
## NkT/V = p_hole/(A_hole*Dt) => p_hole = N*k*T*A_hole*Dt/V
## with v0=0, Dv = (N*k*T*A_hole*Dt)/(V*m) => Dv/Dt = <a> = (N*k*T*A_hole)/(v*m)
## ~ dv/dt = a = (N*k*T*A_hole)/(V*m)
#def accelerate(N,A,V,m):
#    # acelleration per timestep for 1 box
#    return (N*k*T*A)/(V*m)
#m = m_sat
#A_hole = (L/2.)**2
#V = L**3
#a = accelerate(N,A_hole, V, m)
#print 'analytical acceleration per box per timestep'
#print a
#print '---------------------------------------------------------'
#
#
#
#
#	### 2) w/ momentum and mass loss from rocket engine. simulate full acceleration
#		# suitable timstep between 0 and 20 mins (1200?). total mass = m_sat + m_full
#		# m_full is what we calculated last exercise. 
#		# eval. Dv per timestep + updated mass.
#		# what speed after 20 mins? rel error from v_esc? 
#
## number of particles out in 1ns: 8457198. n_o = 8457198. dn/dt = n_o/Dt.
## dm_o/dt = n_o/Dt*m_h
#n_o = 8457198 # number of particles out in 1e-9 seconds. 
#dm = n_o*m_h/Dt
#timespan = 1200 # s timespan to reach v_max
#dt = 1 # s, timestep per interaction
#vel = 0
#m_full = 30387.6201246 # kg. mass of fuel to reach escape vel. in 20 mins w/out fuel mass
#nr_box = 894436114866 # nr. of boxes w/o/ fuel mass.
#m = m_sat+m_full
#for timestep in range(timespan):
#    m += - nr_box*dm*dt
#    if m == m_sat:
#	print 'fuel empty after',timestep,'seconds'
#    a = nr_box*accelerate(N,A_hole,V,m)
#    vel += a*dt
#print 'adjusted velocity with fuel mass'
#print vel
#print 'escape velocity'
#print v_esc
#print 'percentage of v_esc', 1-(abs(vel-v_esc)/v_esc)
#print '---------------------------------------------------------'
#
"""
adjusted velocity with fuel mass
421.670349011
escape velocity
16926.307383
percentage of v_esc 0.0249121287632
"""

	### 3) adjust density (size of box), temperature, amount of fuel and nr. of 
		# boxes to reach v_esc in 20 mins. 
#try L = 1e-7
from Rocket import *
L = 1e-7
N = 1e5
r,v = generate(0, L, N, T, m_h, seed=100)
#test_distribution(r,v)
test_values(v,N,T,m_h)
"""
---------------------------------------------------------
kinetic energy
numerical 2.06135231678e-14
analytical 2.07097278e-14
relative error 0.00464538371288
---------------------------------------------------------
velocity
analytical 1024806365.09
numerical 1022554325.78
relative error 0.00464538371288
---------------------------------------------------------
"""
Dt = 1e-9 #s integration time 1 nanosecond.
n = 10000 # 1e4 timesteps during integration.
l_s, l_e  = L/4., 3.*L/4. # begining and end of hole
#print abs(v)
p, n_p = simulate(Dt, n, r, v, L,m_h,  l_s, l_e)






######################################################################
######################################################################
######################################################################
                ####printing runtime:####
print '---------------------------------------------------------'
print 'runtime'
print time.clock() - start_time, 'seconds'
print '---------------------------------------------------------'
"""
version: 0.1 | 170916
notes: make functions and classes for 1A.6 later, to be callable here. 
"""
