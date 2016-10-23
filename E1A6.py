"""
Exercise 1A.6: design rocket engine for flight of satelite of 
"""
from AST1100SolarSystem import AST1100SolarSystem
from numpy import(
pi, sqrt, random, zeros, linspace, fft, histogram, mean, sum, where, nonzero, shape
)
from matplotlib.pyplot import plot, show, figure, bar, title, savefig
from MySolarSystem import * # module for all SolarSystem functions.
import time

#	#noting processor clock 
start_time = time.clock()

#	#retrieving solar system data
seed = Myseed() # fetching my seed from textfile myseed.txt through module. 
m_0 = p_mass(seed, 0) #mass of home planet (home is 0)
R = p_radius(seed, 0) #radius of home planet (home is 0) in km.

	#setting important constants
AU = 149597870700 #m exactly defined [2]
G = 6.67408e-11 # N*m^2*kg^-2 Newton's gravitational constant [1]
mS = 1988500e24 #kg solar mass[3]
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
R = R*1000
v_esc = sqrt((2*G*m_p)/R) # m/s
print '---------------------------------------------------------'
print 'escape velocity is:', v_esc, 'm/s' 
print '---------------------------------------------------------'
"""
escape velocity is: 16926.307383 m/s
"""

	#### 2: dedicate vectors for particles in box ####
random.seed(100) # seeding the generator to retain numbers between runs
l = 0
h = L
siz = N
r = random.uniform(l,h,(siz,3))

mu = 0
o = sqrt((k*T)/m_h) # sigma, the standard deviation
v = random.normal(mu,o,(siz,3))


	#### test to check uniform and normal distribution ###
	#hist, bins = histogram(r, 'auto')
	#width = 0.7*(bins[1] - bins[0])
	#center = (bins[:-1]+bins[1:])/2
	#bar(center, hist, align='center', width=width)
	#title('position')
	#savefig('testfigures/histogram_test_position.pdf')
	#show()
	#hist, bins = histogram(v, 'auto')
	#width = 0.7*(bins[1] - bins[0])
	#center = (bins[:-1]+bins[1:])/2
	#bar(center, hist, align='center', width=width)
	#title('velocity')
	#savefig('testfigures/histogram_test_velocity.pdf')
	#show()
	#### Test successfull. r is uniform and v is gaussian ###


print 'testing generated values'
print '---------------------------------------------------------'

	#### 3a: is Ek = (3/2)*k*T? ####

Ek_p = (3./2)*N*k*T
Ek = (0.5*m_h)*sum( v**2 )
rel_err = abs( (Ek - Ek_p)/Ek_p )
print 'kinetic energy'
print 'numerical', Ek_p
print 'analytical', Ek
print 'rel_err', rel_err
print '---------------------------------------------------------'

 
	#### 3b:test that the mean abs velocity follows it's relation ####

absv_b = 2*N*sqrt( (2*k*T)/(pi*m_h) ) # <v> = integ (v*P(v)) | scaled up to attempt to reduce roundoff
absv = sum(sqrt( sum( v**2, axis=1 ) )) #total absolute velocity
relerr_absv = abs((absv- absv_b) / absv_b)
print 'velocity'
print 'analytical', absv_b
print 'numerical', absv
print 'rel_error', relerr_absv
print '---------------------------------------------------------'


	### 4: simulate particle motion over ns with elastic collision of walls ###
print '---------------------------------------------------------'
print 'simulating particle motion in time'
print '---------------------------------------------------------'
Dt = 1e-9 #s, 1ns total time elapsed
n = 1000 # number of timesteps within 1ns.
dt = Dt/n # timestep
p_out = 0
n_out = 0
p_z = 0
for timestep in range(n):
	for xyz in range(3): 
		roof = where(r[:,xyz] > L, 1,0)
		floor = where(r[:,xyz] < 0, 1,0)
		crash = roof+floor
		if xyz == 2: 
    		    # placing hole in floor of z
		    n_x = where(
		    (r[nonzero(floor),0]>L/4.) & (r[nonzero(floor),0]<3.*L/4), 1,0)
		    n_y = where(
		    (r[nonzero(floor),1]>L/4.) & (r[nonzero(floor),1]<3.*L/4), 1,0)
		    n_pout = n_x*n_y
		    n_pin = where(n_x*n_y > 0, 0,1)
		    p_out += 2.*m_h*sum(abs(v[nonzero(n_pout),2]))
		    n_out += sum(n_pout)
		    r[nonzero(n_pout), xyz] += L
		    r[nonzero(n_pin),xyz] = - r[nonzero(n_pin),xyz]
		    r[nonzero(roof), xyz] = 2*L - r[nonzero(roof),xyz]
		    v[nonzero(n_pin),xyz] = -v[nonzero(n_pin),xyz]
		else:
		    v[nonzero(crash),xyz] = -v[nonzero(crash),xyz]
		    r[nonzero(roof), xyz] = 2*L - r[nonzero(roof),xyz]
		    r[nonzero(floor), xyz] = -r[nonzero(floor), xyz]
	p_z += 2.*m_h*sum(abs(v[nonzero(roof),2]))
	r += v*dt

print 'n_pout'
print n_pout
print 'nonzero'
print nonzero(n_pout)
print 'v[:,2]'
print v[:,2]


print 'motion simulation'
print '---------------------------------------------------------'
# loop above uses: n, r, L, m_h, v
# and funcs: for, range, where, +, if, ==, nonzero, &, *, >, <, np.sum, abs, else, 

	### 5: choose a wall. count particles that collide over Dt. Count momentum. ###
			# (see above) #

	### a) calculate gas pressure in box w/ momentum and timeperiod ###
		# using pressure P = (F/A) = (1/A)*(dp/dt) => p_z/(A*Dt)
print 'pressure in box'
A = L**2
P_num = p_z/(float(A*Dt)) # Pa, pascal = kg*m^-1*s^-2
print 'numerical', P_num

	### b) using the equation of state ###
		# P = nkT, n is particle density, k is boltzmann and T is temp.
		# to get n, we use our particle number N and divide by the colume L**3
P_ana = N*k*T/(L**3) # Pa, pascal = J*m^-3
print 'analytical', P_ana
rel_err_P = abs(P_num - P_ana) / P_ana
print 'relative error', rel_err_P
print '---------------------------------------------------------'
"""
rel_error 0.00232435501734
numerical 13643.0131419
analytical 13806.4852
relative error 0.0118402370841
With an error of roughly 1%, I assume my simulation is close enough. 
"""

	### 6) Make hole in box. quadratic hole L/2 lengths. count particles
		# count momentum through hole. and number of particles escaping. 
		# Keep box density constant. refill at top when particle leaves. 

print 'momentum out', p_out
print 'particles out', n_out
print '---------------------------------------------------------'
"""
momentum out 4.98691138257e-16
particles out 8457198
"""

	### 7) use momentum loss to calculate speed gain Dv of box(and thereby sattelite) 
		# after period Dt.
		# p_out over 1e-9 seconds. assuming Dv is in m/s, we want dp to be in 
		# s as well, not nanoseconds. Thus, I find dp = p_out/Dt => p_out*1e9*s^-1
		# p/s -> (p/s)*Dt = Dp => Dp/m_sat = Dv.
dp = p_out/Dt
def delta_p(dp, Dt):
	return dp*Dt
del_t = 1200 # 20 min in seconds
Dp = delta_p(dp, del_t)
Dv = Dp/m_sat
print 'speed gain of sattelite per box:', Dv, 'm/s'
print '---------------------------------------------------------'
"""
speed gain of sattelite per box: 5.98429365908e-07 m/s
"""

	### 8) minimum number of boxes to accelerate rocket to escape vel. 
		#in 20 mins?
		#Dv -> (v/box) -> v_esc / Dv ~ n boxes
nr_boxes = int(v_esc/Dv)
print 'boxes required to accelerate sattelite to escape velocity:'
print nr_boxes
print '---------------------------------------------------------'


	### 9) fuel amount necessary to accelerate.
		#we know nr. of particles who left with p_out. know how many left
		#dm/dt = n_out*m_h / Dt => (dm/dt)*20 min = Dm per box
		#Dm per box * nr of boxes = fuel spent
dm = n_out*m_h/Dt
def delta_m(dm, Dt):
	return dm*del_t
Dm_box = delta_m(dm,del_t)
Dm_esc = Dm_box*nr_boxes

print 'fuel consumed for only sattelite to reach v_esc:'
print Dm_esc, 'kg'


######################################################################
		####printing runtime:####
print '---------------------------------------------------------'
print 'runtime'
print time.clock() - start_time, 'seconds'	
print '---------------------------------------------------------'
##############bibliography##########################
"""
#[1] used:31/08-16. Peter J. Mohr, David B. Newell, Barry N. Taylor. 21/07-15. 
#"CODATA Recommended Values of the Fundamental Physical Constants: 2014". 
#rXiv:1507.07956 [physics.atom-ph]
#[2] used:31/08-16. International Astronomical Union, ed. 31/08-12. 
#"RESOLUTION B2 on the re-definition of the astronomical unit of length".
#Beijing, China 
#[3] used:31/08-16. Dr. David R. Williams. 29/02-16. "Sun Fact Sheet". 
#NASA Goddard Space Flight Center. 
#"http://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html"
[4] used 01/09-16. Dr John Huchra. ~2005. "Physical and astrnomical constants".
https://www.cfa.harvard.edu/~dfabricant/huchra/ay145/constants.html
"""
